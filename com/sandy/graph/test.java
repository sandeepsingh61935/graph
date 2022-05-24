package com.sandy.graph;
public class test {
    public static void main(String[] args) {
        Graph g = new Graph(6);
        g.addEdge(0,1);
        g.addEdge(0,2);
        g.addEdge(1,3);
        g.addEdge(2,4);
        g.addEdge(4,5);
        g.addEdge(5,3);
        g.addEdge(1,5);
        g.addEdge(2,5);
        g.addEdge(0,5);
        g.addEdge(0,4);
        g.addEdge(0,5);
        g.addEdge(1,2);
        g.addEdge(2,3);
        g.addEdge(1,3);
        g.addEdge(1,5);
        g.addEdge(2,3);
        g.addEdge(2,5);
        g.addEdge(3,0);
        g.addEdge(3,1);
        g.addEdge(3,2);
        g.addEdge(3,4);
        g.addEdge(4,0);
        g.addEdge(4,1);
        g.addEdge(4,2);
        g.addEdge(4,3);
        g.addEdge(4,5);
        g.addEdge(5,1);
        g.addEdge(5,2);
        g.addEdge(5,3);
        g.addEdge(5,4);
        g.addEdge(5,0);
        g.printGraph();
        System.out.println("BFS of graph is: " + g.bfsGraph());
        System.out.println("DFS of graph is: " + g.dfsGraph());
        System.out.println("Detect cycle in graph: " + g.detectCycleDFS());
        System.out.println("Detect cycle in graph: " + g.detectCycleBFS());
        System.out.println("Detect cycle in graph: " + g.detectCycleDFSBi());
    }
 }
