use crate::ancestors::AncestralSequence;
use std::io;
use std::io::Write;
use std::path::Path;

/// An interval in an ancestor that is covered by a parent node in the tree sequence.
/// The interval is defined by the start and (exclusive) end position of the interval and the index of the
/// parent node.
#[derive(Debug, Clone)]
pub struct TreeSequenceInterval {
    pub(crate) parent: usize,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

impl TreeSequenceInterval {
    pub fn new(parent: usize, start: usize, end: usize) -> Self {
        Self { parent, start, end }
    }
}

/// A node in the tree sequence. The node is defined by the index of the ancestor sequence it
/// represents and a list of intervals that define what parent nodes cover the ancestor sequence.
#[derive(Debug, Clone)]
pub struct TreeSequenceNode {
    // todo hide fields
    pub(crate) ancestor_index: usize,
    pub(crate) node_intervals: Vec<TreeSequenceInterval>,
}

impl TreeSequenceNode {
    pub fn new(ancestor_index: usize, node_intervals: Vec<TreeSequenceInterval>) -> Self {
        TreeSequenceNode {
            ancestor_index,
            node_intervals,
        }
    }

    pub fn empty(ancestor_index: usize) -> Self {
        TreeSequenceNode {
            ancestor_index,
            node_intervals: Vec::new(),
        }
    }

    pub fn tskit_format_node(
        &self,
        ancestor: &AncestralSequence,
        writer: &mut dyn Write,
    ) -> io::Result<()> {
        writer.write_fmt(format_args!(
            "{id}\t{is_sample}\t{time}\n",
            id = self.ancestor_index + 1, // add one to the node index because tskit uses the virtual root node, so we encode the root twice
            is_sample = 0, // todo samples are not supported yet
            time = ancestor.relative_age(),
        ))
    }

    pub fn tskit_format_edges(&self, writer: &mut dyn Write) -> io::Result<()> {
        for interval in &self.node_intervals {
            writer.write_fmt(format_args!(
                "{left}\t{right}\t{parent}\t{child}\n",
                left = interval.start, // todo get the actual genomic position instead of the index into the variant sites
                right = interval.end,  // same here
                // add one to the node index because tskit uses the virtual root node, so we encode the root twice
                parent = interval.parent + 1,
                child = self.ancestor_index + 1,
            ))?;
        }
        Ok(())
    }
}

// todo this should probably be a proper struct with hidden fields
pub struct TreeSequence(pub Vec<TreeSequenceNode>, pub Vec<AncestralSequence>);

impl TreeSequence {
    pub fn tskit_export(&self, path: &Path) -> io::Result<()> {
        let mut node_file = path.to_path_buf();
        node_file.push("nodes.tsv");
        let mut writer = std::fs::File::create(node_file)?;

        writer.write_fmt(format_args!("id\tis_sample\ttime\n"))?;
        // write root node twice, because tskit uses the virtual root
        writer.write_fmt(format_args!("0\t0\t1.2\n"))?;

        for node in &self.0 {
            node.tskit_format_node(&self.1[node.ancestor_index], &mut writer)?;
        }

        let mut edge_file = path.to_path_buf();
        edge_file.push("edges.tsv");
        let mut writer = std::fs::File::create(edge_file)?;

        writer.write_fmt(format_args!("left\tright\tparent\tchild\n"))?;
        // write edge from virtual root to root
        writer.write_fmt(format_args!("0\t{}\t0\t1\n", self.0[1].node_intervals[0].end))?;

        // skip first edge because tskit doesn't like the root node to have an edge to itself. TODO we can remove this anyway at some point
        for node in self.0.iter().skip(1) {
            node.tskit_format_edges(&mut writer)?;
        }
        Ok(())
    }
}
