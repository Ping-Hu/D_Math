# discrete_math
Homework 02 solution

– Name, ID, email address

– Explain how to use your solution code
The command for running the solution is as follows:
cse547.exe mesh.m

– Source codes
– Project files
– Explain your algorithm for each requirement
(1). Trace Boundaries
Method one: implement boundary tracing by myself. First, find an edge with only one half edge which is a boundary edge. Second, trace this boundary edge to find the adjacent boundary edge by checking if he->he_next()->he_sym() == NULL. Third, store found boundary edges in a list called boundary.
Method two: directly call the method bound.loops().size()
(2). Euler Number
Euler number is computed according to the function (numFaces + numVertices -numEdges).
(3). Face Normal
d = (v1 → point() − v0 → point) × (v2 → point() − v0 → point)

(4). Vertex Normal

(5). Vertex Gaussian Curvature

(6). Gauss-Bonnet Formula


Comments:
1. The only edited files are main.cpp and CurvatureMesh.h.
2. Build a new project with the same relative path of "MeshLib", "CurvatureMesh" and the file, main.cpp.
3. I don't think my code is very compact and efficient. Specifically, there must be some easier ways to trace boundaries. But I didn't go through all Gu's codes for which reason I only use the most basic part of the data structure in his codes when implementing boundary tracing.
