package com.sandy.graph;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Queue;

class Pair<T1, T2> {
    private T1 first;
    private T2 second;

    public Pair(T1 first, T2 second) {
        this.first = first;
        this.second = second;
    }

    public T1 getFirst() {
        return first;
    }

    public T2 getSecond() {
        return second;
    }

    public void setFirst(T1 first) {
        this.first = first;
    }

    public void setSecond(T2 second) {
        this.second = second;
    }

    public String toString() {
        return "(" + first + ", " + second + ")";
    }
}
public class Graph {
    private final int NUM_OF_VERTICES;
    private ArrayList<ArrayList<Integer>> graph;
    private ArrayList<ArrayList<Pair<Integer, Integer>>> weightedGraph;
    public Graph(int numVertices,boolean isWeightedGraph) {
        this.NUM_OF_VERTICES = numVertices;
        if (isWeightedGraph) {
            this.weightedGraph = new ArrayList<ArrayList<Pair<Integer, Integer>>>(this.NUM_OF_VERTICES);
            for (int i = 0; i < this.NUM_OF_VERTICES; i++) {
                this.weightedGraph.add(new ArrayList<>());
            }
        } else {
            this.graph = new ArrayList<ArrayList<Integer>>(this.NUM_OF_VERTICES);
            for (int i = 0; i < this.NUM_OF_VERTICES; i++) {
                this.graph.add(new ArrayList<>());
            }
        }
        
    }
    //add an edge from vertex u to vertex v
    public void addEdge(int u, int v) {
        graph.get(u).add(v);
    }

    // add weighted edge from vertex u to vertex v
    public void addWeightedEdge(int u, int v, int weight) {
        Pair<Integer, Integer> pair = new Pair<>(v, weight);
        weightedGraph.get(u).add(pair);
    }

    // print adjacency list representation of graph
    public void printGraph() {
        for (int i = 0; i < this.NUM_OF_VERTICES; i++) {
            System.out.print(i + " -> ");
            for (int j = 0; j < this.graph.get(i).size(); j++) {
                System.out.print(this.graph.get(i).get(j) + " ");
            }
            System.out.println();
        }
    }

    //return the list of vertices adjacent to vertex v
    public ArrayList<Integer> getAdjacentVertices(int v) {
        return graph.get(v);
    }

    // return the list of weights adjacent to vertex v
    public ArrayList<Pair<Integer, Integer>> getAdjacentVerticesWithWeights(int v) {
        return weightedGraph.get(v);
    }

    //return the number of vertices in the graph
    public int getNumOfVertices() {
        return NUM_OF_VERTICES;
    }

    /**
     * @description breadthfirst search of vertices in the graph
     */
    public ArrayList<Integer> bfsGraph() {
        ArrayList<Boolean> visited = new ArrayList<>();
        ArrayList<Integer> queue = new ArrayList<>();
        ArrayList<Integer> bfs = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            visited.add(false);
        }

        visited.set(0,true);
        queue.add(0);

        while(!queue.isEmpty()) {
            int u = queue.remove(0);
            bfs.add(u);
            ArrayList<Integer> adjacentVertices = getAdjacentVertices(u);
            for(int v : adjacentVertices) {
                if(!visited.get(v)) {
                    visited.set(v,true);
                    queue.add(v);
                }
            }
        }
        return bfs;
    }


    /**
     * @description depthfirst search of vertices in the graph
     */
    public ArrayList<Integer> dfsGraph() {
        ArrayList<Boolean> visited = new ArrayList<>();
        ArrayList<Integer> dfs = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            visited.add(false);
        }

        dfsGraphUtil(0, visited, dfs);
        return dfs;
    }

    /**
     * @param i the index of the vertex to be visited
     * @param dfs the dfs list
     * @param visited the list of vertices visited
     * @description recursive function to perform depthfirst search
     */
    private void dfsGraphUtil(int i, ArrayList<Boolean> visited, ArrayList<Integer> dfs) {
        visited.set(i,true);
        dfs.add(i);
        ArrayList<Integer> adjacentVertices = getAdjacentVertices(i);
        for(int v : adjacentVertices) {
            if(!visited.get(v)) {
                dfsGraphUtil(v, visited, dfs);
            }
        }
    }


    // detect cycle in bidirected graph
    public boolean detectCycleDFSBi() {
        ArrayList<Boolean> visited = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            visited.add(false);
        }

        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            if (!visited.get(i)) {
                if (detectCycleDFSUtilBi(i, visited)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * @param i the index of the vertex to be visited
     * @param visited the list of vertices visited
     * @description recursive function to perform depthfirst search
     */
    private boolean detectCycleDFSUtilBi(int i, ArrayList<Boolean> visited) {
        visited.set(i,true);
        ArrayList<Integer> adjacentVertices = getAdjacentVertices(i);
        for(int v : adjacentVertices) {
            if(!visited.get(v)) {
                if(detectCycleDFSUtilBi(v, visited)) {
                    return true;
                }
            } else {
                return true;
            }
        }
        visited.set(i,false);
        return false;
    }

    // detect cycle in directed graph using BFS algorithm
    // will use KAHN's algorithm which only supports directed acyclic graph
    //  decising condition: topSort length == NUM_OF_VERTICES return false otherwise return true
    public boolean detectCycleUsingBFS() {
        // indegree of each vertex is stored in this list  
        ArrayList<Integer> indegree = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            for(int  v : getAdjacentVertices(i)) {
                indegree.set(v, indegree.get(v) + 1);
            }
        }
        Queue<Integer> q = new LinkedList<>();
        // vertices with indegree 0 are added to the queue becoz they are the starting vertices(top of the graph)
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(indegree.get(i) == 0) {
                q.add(i);
            }
        }
        int cnt = 0;
        while(!q.isEmpty()) {
            int u = q.remove();
            cnt++;
            for(int v : getAdjacentVertices(u)) {
                indegree.set(v, indegree.get(v) - 1);
                if(indegree.get(v) == 0) {
                    q.add(v);
                }
            }
        }
        return cnt != NUM_OF_VERTICES ;
    }
    // detectCycle using BFS algorithm
    public boolean detectCycleBFS()  {
        ArrayList<Boolean> visited = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            visited.add(false);
        }

        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            if (!visited.get(i)) {
                if (detectCycleBFSUtil(i, visited)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * @description util method to detect cycle using BFS
     * @param i the index of the vertex to be visited
     * @param visited the list of vertices visited
     * @return  true if cycle is detected, false otherwise
     */
    private boolean detectCycleBFSUtil(int i, ArrayList<Boolean> visited) {
        ArrayList<Pair<Integer, Integer>> queue = new ArrayList<>();
        queue.add(new Pair<>(i,-1));
        visited.add( i, true);
        while(!queue.isEmpty()) {
            Pair<Integer, Integer> pair = queue.remove(0);
            int u = pair.getFirst();
            int parent = pair.getSecond();
            ArrayList<Integer> adjacentVertices = getAdjacentVertices(u);
            for(int v : adjacentVertices) {
                if(!visited.get(v)) {
                    visited.set(v,true);
                    queue.add(new Pair<>(v,u));
                } else if(v != parent) {
                    return true;
                }
            }
        }
        return false;
    }

    // detect cycle using DFS in undiredected graph 
    public boolean detectCycleDFS(){
        ArrayList<Boolean> visited = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            visited.add(false);
        }

        for(int i =0 ;i < visited.size(); i++) {
            if(!visited.get(i)) {
                if(detectCycleDFSUtil(i,-1, visited)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * @param i the index of the vertex to be visited
     * @param parent the parent of the vertex
     * @param visited the list of vertices visited
     * @return  true if cycle is detected, false otherwise
     */
    private boolean detectCycleDFSUtil(int i ,int parent, ArrayList<Boolean> visited) {
        visited.set(i,true);
        ArrayList<Integer> adjacentVertices = getAdjacentVertices(i);
        for(int v : adjacentVertices) {
            if(!visited.get(v)) {
                if(detectCycleDFSUtil(v,i, visited)) {
                    return true;
                }
            } else if(v != parent) {
                return true;
            }
        }
        return false;
    }

    // detect graph is bipartite graph or not
    public boolean isBipartiteUsingDFS() {
        ArrayList<Integer> colors =  new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            colors.add(-1);
        }

        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(colors.get(i) == -1) {
                if(!isBipartiteUtilUsingDFS(i, colors)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @param i the index of the vertex to be visited
     * @param colors the list of colors of the vertices
     * @return  true if bipartite graph is detected, false otherwise
     */
    public boolean isBipartiteUtilUsingDFS(int i, ArrayList<Integer> colors) {
        if ( colors.get(i) == - 1) 
            colors.set(i,0);
        ArrayList<Integer> adjacentVertices = getAdjacentVertices(i);
        for(int v : adjacentVertices) {
            if(colors.get(v) == colors.get(i)) {
                return false;
            } else if(colors.get(v) != colors.get(i)) {
                colors.set(v, 1 - colors.get(i));
                if(!isBipartiteUtilUsingDFS(v, colors)) {
                    return false;
                }
            }
        }
        return true;
    }


    // detect graph is bipartite graph or not using BFS
    public boolean isBipartiteGraphUsingBFS() {
        ArrayList<Integer> colors = new ArrayList<Integer>();

        for(int i =0 ;i < NUM_OF_VERTICES;i++) {
            colors.add(-1);
        }

        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(colors.get(i) == -1) {
                if(!isBipartiteUtilUsingBFS(i, colors)) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * @param i the index of the vertex to be visited
     * @param colors the list of colors of the vertices
     * @return  true if the vertex is visited, false otherwise
     */
    private boolean isBipartiteUtilUsingBFS(int i, ArrayList<Integer> colors) {
        ArrayList<Integer> queue = new ArrayList<Integer>();
        queue.add(i);
        colors.set(i,0);

        while(!queue.isEmpty()) {
            int u = queue.remove(0);
            ArrayList<Integer> adjacentVertices = getAdjacentVertices(u);
            for(int v : adjacentVertices) {
                if (colors.get(v) == - 1) {
                    colors.set(v,1 - colors.get(i));
                    queue.add(v);
                }
                else if (colors.get(v) == colors.get(i)) {
                    return false;
                }
            }
        }
        return true;
    }

    // topological sort using DFS
    public ArrayList<Integer> topologicalSortUsingDFS() {
        ArrayList<Boolean> visited = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            visited.add(false);
        }

        ArrayList<Integer> stack = new ArrayList<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(!visited.get(i)) {
                topologicalSortUtilUsingDFS(i, visited, stack);
            }
        }
        Collections.reverse(stack);
        return stack;
    }

    /**
     * @param i the index of the vertex to be visited
     * @param visited the list of vertices visited
     * @param stack the stack of vertices
     */
    private void topologicalSortUtilUsingDFS(int i, ArrayList<Boolean> visited, ArrayList<Integer> stack) {
        visited.set(i,true);
        ArrayList<Integer> adjacentVertices = getAdjacentVertices(i);
        for(int v : adjacentVertices) {
            if(!visited.get(v)) {
                topologicalSortUtilUsingDFS(v, visited, stack);
            }
        }
        stack.add(i);
        return;
    }

    // topological sort using BFS(KAHN's Algorithm)
    // works only for DAG(Directed Acyclic Graph)
    // Space Complexity: O(V+E)
    // Time complexity: O(V+E)
    public ArrayList<Integer> topologicalSortUsingBFS() {
        // indegree of each vertex is stored in this list  
        ArrayList<Integer> indegree = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            for(int  v : getAdjacentVertices(i)) {
                indegree.set(v, indegree.get(v) + 1);
            }
        }
        Queue<Integer> q = new LinkedList<>();
        ArrayList<Integer> topSort = new ArrayList<>();
        // vertices with indegree 0 are added to the queue becoz they are the starting vertices(top of the graph)
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(indegree.get(i) == 0) {
                q.add(i);
            }
        }

        while(!q.isEmpty()) {
            int u = q.remove();
            topSort.add(u);
            for(int v : getAdjacentVertices(u)) {
                indegree.set(v, indegree.get(v) - 1);
                if(indegree.get(v) == 0) {
                    q.add(v);
                }
            }
        }
        return topSort;
    }


    // shortest distance between source and all the vertices USING BFS algorithm
    // ASSUMPTION: unit weights
    // Time Complexity: O(V+E)
    // Space Complexity: O(V)
    public ArrayList<Integer> shortestDistanceUsingBFS(int source) {
        ArrayList<Integer> distance = new ArrayList<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            distance.add(Integer.MAX_VALUE);
        }
        distance.set(source, 0);
        ArrayList<Integer> queue = new ArrayList<>();
        queue.add(source);
        while(!queue.isEmpty()) {
            int u = queue.remove(0);
            ArrayList<Integer> adjacentVertices = getAdjacentVertices(u);
            for(int v : adjacentVertices) {
                if(distance.get(v) == Integer.MAX_VALUE) {
                    distance.set(v, distance.get(u) + 1);
                    queue.add(v);
                }
            }
        }
        return distance;
    }
    
    // shortest distance between source and all the vertices USING Dijkstra's algorithm
    // Time Complexity: O((N+E)logN) = O((NlogN)
    // Space Complexity: O(NUM_OF_VERTICES)
    public ArrayList<Integer> shortestDistanceUsingDijkstra(int source) {
        ArrayList<Integer> distance = new ArrayList<>();
        ArrayList<Boolean> visited = new ArrayList<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            distance.add(Integer.MAX_VALUE);
            visited.add(false);
        }
        distance.set(source, 0);
        PriorityQueue<Pair<Integer,Integer>> pq = new PriorityQueue<>();
        pq.add(new Pair<>(0,source));
        while(!pq.isEmpty()) {
            int v = pq.remove().getSecond();
            visited.set(v, true);
            ArrayList<Pair<Integer,Integer>> adjacentVertices = getAdjacentVerticesWithWeights(v);
            for(Pair<Integer,Integer> pair : adjacentVertices) {
                int u = pair.getSecond();
                int w = pair.getFirst();
                if (!visited.get(u)) {
                    if (distance.get(u) > distance.get(v) + w) {
                        distance.set(u, distance.get(v) + w);
                        pq.add(new Pair<>(distance.get(u), u));
                    }
                }
            }
        }
        return distance;
    }


    

}