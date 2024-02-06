use crate::ancestors::{Ancestor, AncestorArray, AncestralSequence, VariantIndex};
use crate::dna::SequencePosition;
use std::io;
use std::io::Write;
use std::path::Path;

/// An interval in an ancestor that is covered by a parent node in the tree sequence.
/// The interval is defined by the start and (exclusive) end position of the interval and the index of the
/// parent node.
#[derive(Debug, Clone)]
pub struct TreeSequenceEdge {
    pub parent: usize,
    pub start: usize,
    pub end: usize,
}

impl TreeSequenceEdge {
    pub fn new(parent: usize, start: SequencePosition, end: SequencePosition) -> Self {
        Self {
            parent,
            start: start.unwrap(),
            end: end.unwrap(),
        }
    }
}

/// A node in the tree sequence. The node is defined by the index of the ancestor sequence it
/// represents and a list of intervals that define what parent nodes cover the ancestor sequence.
#[derive(Debug, Clone)]
pub struct TreeSequenceNode {
    ancestor_index: usize,
    edges: Vec<TreeSequenceEdge>,
    mutations: Vec<usize>,
}

impl TreeSequenceNode {
    pub(crate) fn new(
        ancestor_index: usize,
        edges: Vec<TreeSequenceEdge>,
        mutations: &[VariantIndex],
    ) -> Self {
        TreeSequenceNode {
            ancestor_index,
            edges,
            mutations: mutations.iter().map(|v| v.unwrap()).collect(),
        }
    }

    pub(crate) fn empty(ancestor_index: usize) -> Self {
        TreeSequenceNode {
            ancestor_index,
            edges: Vec::new(),
            mutations: Vec::new(),
        }
    }

    /// Get the index of the ancestor sequence that this node represents
    pub fn ancestor(&self) -> usize {
        self.ancestor_index
    }

    /// Get the collection of [`TreeSequenceEdges`] that define the intervals covered by parent nodes.
    ///
    /// [`TreeSequenceEdges`]: TreeSequenceEdge
    pub fn edges(&self) -> &[TreeSequenceEdge] {
        &self.edges
    }

    /// Get a collection of mutation sites. The collection is a list of indices into the variant
    /// site vector.
    pub fn mutations(&self) -> &[usize] {
        &self.mutations
    }

    pub fn tskit_format_node(
        &self,
        ancestor: &AncestralSequence,
        writer: &mut dyn Write,
    ) -> io::Result<()> {
        writer.write_fmt(format_args!(
            "{id}\t{is_sample}\t{time}\n",
            id = self.ancestor_index + 1, // add one to the node index because tskit uses the virtual root node, so we encode the root twice
            is_sample = 1,                // todo samples are not supported yet
            time = ancestor.relative_age(),
        ))
    }

    pub fn tskit_format_edges(&self, writer: &mut dyn Write) -> io::Result<()> {
        for interval in &self.edges {
            writer.write_fmt(format_args!(
                "{left}\t{right}\t{parent}\t{child}\n",
                left = interval.start,
                right = interval.end,
                // add one to the node index because tskit uses the virtual root node, so we encode the root twice
                parent = interval.parent + 1,
                child = self.ancestor_index + 1,
            ))?;
        }
        Ok(())
    }

    pub fn tskit_format_mutations(&self, writer: &mut dyn Write) -> io::Result<()> {
        for mutation in &self.mutations {
            writer.write_fmt(format_args!(
                "{site}\t{node}\t{derived_state}\n",
                site = mutation, // add one to the node index because tskit uses the virtual root node, so we encode the root twice
                node = self.ancestor_index + 1,
                derived_state = 'A', // todo get the actual derived state
            ))?;
        }
        Ok(())
    }
}

/// A tree sequence consisting of a set of ancestors and corresponding nodes in the tree sequence.
/// The node indices correspond to the ancestor indices.
pub struct TreeSequence {
    pub nodes: Vec<TreeSequenceNode>,
    pub ancestors: AncestorArray,
}

impl TreeSequence {
    pub fn tskit_export(&self, path: &Path) -> io::Result<()> {
        let mut node_file = path.to_path_buf();
        node_file.push("nodes.tsv");
        let mut writer = std::fs::File::create(node_file)?;

        writer.write_fmt(format_args!("id\tis_sample\ttime\n"))?;
        // write root node twice, because tskit uses the virtual root
        writer.write_fmt(format_args!(
            "{id}\t{is_sample}\t{time}\n",
            id = 0,
            is_sample = 1,
            time = 2.0,
        ))?;

        for node in &self.nodes {
            node.tskit_format_node(&self.ancestors[Ancestor(node.ancestor_index)], &mut writer)?;
        }

        let mut edge_file = path.to_path_buf();
        edge_file.push("edges.tsv");
        let mut writer = std::fs::File::create(edge_file)?;

        writer.write_fmt(format_args!("left\tright\tparent\tchild\n"))?;
        // write edge from virtual root to root
        writer.write_fmt(format_args!("0\t{}\t0\t1\n", self.nodes[1].edges[0].end))?;

        // skip first edge because tskit doesn't like the root node to have an edge to itself. TODO we can remove this anyway at some point
        for node in self.nodes.iter().skip(1) {
            node.tskit_format_edges(&mut writer)?;
        }

        let mut mutation_file = path.to_path_buf();
        mutation_file.push("mutations.tsv");
        let mut writer = std::fs::File::create(mutation_file)?;

        writer.write_fmt(format_args!("site\tnode\tderived_state\n"))?;
        for node in &self.nodes {
            node.tskit_format_mutations(&mut writer)?;
        }
        Ok(())
    }
}
