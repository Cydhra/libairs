use crate::ancestors::{Ancestor, AncestorArray, AncestralSequence};
use crate::variants::SequencePosition;
use std::io;
use std::io::{BufWriter, Write};
use std::path::Path;

/// A mutation in an ancestor sequence. The mutation is defined by the sequence position
/// and the derived state in FASTA notation.
#[derive(Debug, Clone, Hash, Eq, PartialEq)]
pub struct Mutation {
    pub variant_index: u32,
    pub derived_state: char,
}

/// An interval in an ancestor that is covered by a parent node in the tree sequence.
/// The interval is defined by the start and (exclusive) end position of the interval and the index of the
/// parent node.
#[derive(Debug, Clone, PartialEq)]
pub struct TreeSequenceEdge {
    pub parent: u32,
    pub start: u32,
    pub end: u32,
}

impl TreeSequenceEdge {
    pub fn new(parent: u32, start: SequencePosition, end: SequencePosition) -> Self {
        Self {
            parent,
            start: start.unwrap(),
            end: end.unwrap(),
        }
    }
}

/// A node in the tree sequence. The node is defined by the index of the ancestor sequence it
/// represents and a list of intervals that define what parent nodes cover the ancestor sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct TreeSequenceNode {
    ancestor_index: u32,
    edges: Vec<TreeSequenceEdge>,
    mutations: Vec<Mutation>,
    is_sample: bool,
}

impl TreeSequenceNode {
    pub(crate) fn new(
        ancestor_index: u32,
        edges: Vec<TreeSequenceEdge>,
        mutations: Vec<Mutation>,
        is_sample: bool,
    ) -> Self {
        TreeSequenceNode {
            ancestor_index,
            edges,
            mutations,
            is_sample,
        }
    }

    /// Get the index of the ancestor sequence that this node represents
    pub fn ancestor(&self) -> usize {
        self.ancestor_index as usize
    }

    /// Get the collection of [`TreeSequenceEdges`] that define the intervals covered by parent nodes.
    ///
    /// [`TreeSequenceEdges`]: TreeSequenceEdge
    pub fn edges(&self) -> &[TreeSequenceEdge] {
        &self.edges
    }

    /// Get a collection of mutation sites. The collection is a list of indices into the variant
    /// site vector.
    pub fn mutations(&self) -> &[Mutation] {
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
            is_sample = self.is_sample as u8,
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
                site = mutation.variant_index, // add one to the node index because tskit uses the virtual root node, so we encode the root twice
                node = self.ancestor_index + 1,
                derived_state = mutation.derived_state,
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
        let mut writer = BufWriter::new(std::fs::File::create(node_file)?);

        writer.write_fmt(format_args!("id\tis_sample\ttime\n"))?;
        // write root node twice, because tskit uses the virtual root
        writer.write_fmt(format_args!(
            "{id}\t{is_sample}\t{time}\n",
            id = 0,
            is_sample = 0,
            time = 3.0,
        ))?;

        for node in &self.nodes {
            let sample_data = self.ancestors.generate_sample_data();
            if !node.is_sample {
                node.tskit_format_node(
                    &self.ancestors[Ancestor(node.ancestor_index)],
                    &mut writer,
                )?;
            } else {
                node.tskit_format_node(
                    &sample_data[node.ancestor_index as usize - self.ancestors.len() as usize],
                    &mut writer,
                )?;
            }
        }

        let mut edge_file = path.to_path_buf();
        edge_file.push("edges.tsv");
        let mut writer = BufWriter::new(std::fs::File::create(edge_file)?);

        writer.write_fmt(format_args!("left\tright\tparent\tchild\n"))?;
        // write edge from virtual root to root
        writer.write_fmt(format_args!("0\t{}\t0\t1\n", self.nodes[1].edges[0].end))?;

        // skip first edge because tskit doesn't like the root node to have an edge to itself. TODO we can remove this anyway at some point
        for node in self.nodes.iter().skip(1) {
            node.tskit_format_edges(&mut writer)?;
        }

        let mut mutation_file = path.to_path_buf();
        mutation_file.push("mutations.tsv");
        let mut writer = BufWriter::new(std::fs::File::create(mutation_file)?);

        writer.write_fmt(format_args!("site\tnode\tderived_state\n"))?;
        for node in &self.nodes {
            node.tskit_format_mutations(&mut writer)?;
        }
        Ok(())
    }
}
