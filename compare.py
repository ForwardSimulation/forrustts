import tskit

a = tskit.load("regular.trees")
b = tskit.load("consume.trees")

print(a.num_nodes, b.num_nodes)
print(a.num_edges, b.num_edges)

assert a == b
