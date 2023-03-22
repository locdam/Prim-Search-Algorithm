# Prim-s-Search
Generate a minimum spanning tree using Prim's Algorithm

A minimum spanning tree T(V’, E’) is a subset of graph G(V, E) with the same number of vertices as of graph G (V’ = V) and edges equal to the number of vertices of graph G minus one (E’ = |V| - 1). Prim's approach identifies the subset of edges that includes every vertex in the graph, and allows the sum of the edge weights to be minimized.

Prim’s algorithm starts with a single node and works its way through several adjacent nodes, exploring all of the connected edges along the way. Edges with the minimum weights that do not cause cycles in the graph get selected for t inclusion in the MST structure. 

Let’s first look into the steps involved in Prim’s Algorithm to generate a minimum spanning tree:

Step 1: Determine the arbitrary starting vertex.
Step 2: Keep repeating steps 3 and 4 until the fringe vertices (vertices not included in MST)remain. 
Step 3: Select an edge connecting the tree vertex and fringe vertex having the minimum weight.
Step 4: Add the chosen edge to MST if it doesn’t form any closed cycle.
Step 5: Exit
