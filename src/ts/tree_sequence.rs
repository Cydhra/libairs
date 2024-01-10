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

    pub fn tskit_format_node(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_fmt(format_args!(
            "{id}\t{is_sample}\t{time}\n",
            id = self.ancestor_index,
            is_sample = 0, // todo samples are not supported yet
            time = 0.0, // todo get the relative time of the ancestor
        ))
    }

    pub fn tskit_format_edges(&self, writer: &mut dyn Write) -> io::Result<()> {
        for interval in &self.node_intervals {
            writer.write_fmt(format_args!(
                "{left}\t{right}\t{parent}\t{child}\n",
                left = interval.start, // todo get the actual genomic position instead of the index into the variant sites
                right = interval.end, // same here
                parent = interval.parent,
                child = self.ancestor_index,
            ))?;
        }
        Ok(())
    }
}

// todo this should probably be a proper struct with hidden fields
pub struct TreeSequence(pub Vec<TreeSequenceNode>);

impl TreeSequence {
    pub fn tskit_export(&self, path: &Path) -> io::Result<()> {
        let mut node_file = path.to_path_buf();
        node_file.push("nodes.tsv");
        let mut writer = std::fs::File::create(node_file)?;

        writer.write_fmt(format_args!(
            "id\tis_sample\ttime\n"
        ))?;
        for node in &self.0 {
            node.tskit_format_node(&mut writer)?;
        }

        let mut edge_file = path.to_path_buf();
        edge_file.push("edges.tsv");
        let mut writer = std::fs::File::create(edge_file)?;

        writer.write_fmt(format_args!(
            "left\tright\tparent\tchild\n"
        ))?;
        for node in &self.0 {
            node.tskit_format_edges(&mut writer)?;
        }
        Ok(())
    }
}