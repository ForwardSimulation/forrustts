#[cfg(test)]
mod test {
    use crate::wright_fisher::*;

    //FIXME: remove these when done
    //use crate::simplify_tables::simplify_tables;
    //use std::fs::File;
    //use std::io::Write;

    #[test]
    fn test_wright_fisher() {
        let mut tables = neutral_wf(42, 250, 5000, 1000000, 5e-3, 0.0, Some(100));

        //let mut nodes_file = File::create("nodes.txt").unwrap();
        //for n in tables.nodes().iter() {
        //    write!(nodes_file, "{}\n", -1 * (n.time - 5000)).unwrap();
        //}
        //let mut edges_file = File::create("edges.txt").unwrap();
        //for e in tables.edges().iter() {
        //    write!(
        //        edges_file,
        //        "{} {} {} {}\n",
        //        e.left, e.right, e.parent, e.child
        //    )
        //    .unwrap();
        //}

        //let mut samples = Vec::<i32>::new();
        //for (i, n) in tables.nodes().iter().enumerate() {
        //    if n.time == 5000 {
        //        samples.push(i as i32);
        //    }
        //}
        //tables.sort_tables_for_simplification();
        //simplify_tables(&samples, &mut tables);
        //let mut nodes_file = File::create("nodes_simplified.txt").unwrap();
        //for n in tables.nodes().iter() {
        //    write!(nodes_file, "{}\n", -1 * (n.time - 5000)).unwrap();
        //}
        //let mut edges_file = File::create("edges_simplified.txt").unwrap();
        //for e in tables.edges().iter() {
        //    write!(
        //        edges_file,
        //        "{} {} {} {}\n",
        //        e.left, e.right, e.parent, e.child
        //    )
        //    .unwrap();
        //}
    }
}
