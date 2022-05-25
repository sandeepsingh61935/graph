package com.sandy.graph;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Queue;

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
    // (u--> [(v,weight)])
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

    // transpose the graph
    public Graph transposeGraph(boolean isWeightedGraph) {
        Graph g = new Graph(this.NUM_OF_VERTICES, isWeightedGraph);
        for (int i = 0; i < this.NUM_OF_VERTICES; i++) {
            for (int j = 0; j < this.graph.get(i).size(); j++) {
                g.addEdge(this.graph.get(i).get(j), i);
            }
        }
        return g;
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
    //topological sort using DFS in weighted DAG
    public ArrayList<Integer> weightedTopologicalSortUsingDFS() {
        ArrayList<Boolean> visited = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            visited.add(false);
        }

        ArrayList<Integer> stack = new ArrayList<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(!visited.get(i)) {
                weightedTopologicalSortUtilUsingDFS(i, visited, stack);
            }
        }
        return stack;
    }

    // weightedTopologicalSortUtilUsingDFS 
    private void weightedTopologicalSortUtilUsingDFS(int i, ArrayList<Boolean> visited, ArrayList<Integer> stack) {
        visited.set(i,true);
        ArrayList<Pair<Integer, Integer>> adjacentVertices = this.getAdjacentVerticesWithWeights(i);
        for(Pair<Integer, Integer> pair : adjacentVertices) {
            if(!visited.get(pair.getFirst())) {
                weightedTopologicalSortUtilUsingDFS(pair.getFirst(), visited, stack);
            }
        }
        stack.add(i);
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
    /**
     * @read
     * shortest distance between source and all the vertices
        graph type:DAG
        Time Complexity: O(V+E)
        Space Complexity: O(V) 
        steps: 
            1. topSort
            2. normal bfs traversal
            3. update distance
        @param source the source vertex
     */

    public ArrayList<Integer> shortestDistanceUsingBFSForDAG(int source) {
        ArrayList<Integer> distance = new ArrayList<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            distance.add(Integer.MAX_VALUE);
        }
        distance.set(source, 0);
        ArrayList<Integer> stack = this.weightedTopologicalSortUsingDFS();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            int u = stack.get(i);
            ArrayList<Pair<Integer,Integer>> adjacentVertices = this.getAdjacentVerticesWithWeights(u);
            for(Pair<Integer,Integer> pair : adjacentVertices) {
                int v = pair.getFirst();
                int w = pair.getSecond();
                if(distance.get(v) == Integer.MAX_VALUE) {
                    distance.set(v, distance.get(u) + w);
                }
                else if(distance.get(v) > distance.get(u) + w) {
                    distance.set(v, distance.get(u) + w);
                }
            }
        }
        return distance;
    }

    // TC : O(VLOGV + E)
    // SC : O(V)
    public void minimumSpanningTree() {
        ArrayList<Integer> distance = new ArrayList<Integer>();
        ArrayList<Boolean> isMST = new ArrayList<> ();
        ArrayList<Integer> parent = new ArrayList<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            distance.add(Integer.MAX_VALUE);
            isMST.add(false);
            parent.add(-1);
        }

        distance.set(0,0);
        parent.set(0,-1);
        // pair<distance, vertex>
        PriorityQueue<Pair<Integer,Integer>> pq = new PriorityQueue<>();
        pq.add(new Pair<>(distance.get(0),0));
        while(!pq.isEmpty()) {
            Pair<Integer,Integer> pair = pq.remove();
            int u = pair.getSecond();
            isMST.set(u, true);
            ArrayList<Pair<Integer,Integer>> adjacentVertices = this.getAdjacentVerticesWithWeights(u);
            for(Pair<Integer,Integer> pair1 : adjacentVertices) {
                int v = pair1.getFirst();
                int w = pair1.getSecond();
                if(!isMST.get(v) && distance.get(v) > w) {
                    parent.set(u,v);
                    distance.set(v, w);
                    pq.add(new Pair<>(w, v));
                }
            }
        }
        // print spanning tree  
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(parent.get(i) != -1) {
                System.out.println(i + " " + parent.get(i));
            }
        }
    }


    // kruskal's algorithm
    // TC : O(ElogE)
    // SC : O(E)
    // implement the algorithm
    public void kruskal() {
        PriorityQueue<Pair<Integer,Pair<Integer,Integer>>> pq = new PriorityQueue<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            ArrayList<Pair<Integer,Integer>> adjacentVertices = this.getAdjacentVerticesWithWeights(i);
            for(Pair<Integer,Integer> pair : adjacentVertices) {
                int v = pair.getFirst();
                int w = pair.getSecond();
                pq.add(new Pair<>(w, new Pair<>(i, v)));
            }
        }
        int i = 0;
        int[] parent = new int[NUM_OF_VERTICES];
        for(int j = 0; j < NUM_OF_VERTICES; j++) {
            parent[j] = -1;
        }
        while(!pq.isEmpty() && i < NUM_OF_VERTICES - 1) {
            Pair<Integer,Pair<Integer,Integer>> pair = pq.remove();
            int u = pair.getSecond().getFirst();
            int v = pair.getSecond().getSecond();
            if(parent[u] == -1 && parent[v] == -1) {
                parent[u] = v;
                parent[v] = u;
                i++;
            }
        }
        // print spanning tree  
        for(int j = 0; j < NUM_OF_VERTICES; j++) {
            if(parent[j] != -1) {
                System.out.println(j + " " + parent[j]);
            }
        }
    }


    /**
     * @description implement krushal's algorithm (disjoint set)
     * steps:
     * 1. sort the edges in non-decreasing order of their weight
     * 2. create a disjoint set
     * 3. for each edge, if the two vertices are not in the same set, add the edge to the MST
     * 
     * FINALLY NOTE:
     *  TC: O(ElogE) for sorting the edges + O(E*4C) for creating the disjoint set C is the constant values.so, O(ELOGE)    
     *  SC: O(E) sorting , O(V) for parent and O(V) for rank.O(V+E)
     */

     public void MSTUsingKruskals() {
         // sort the edges in non-decreasing order of their weight
         // <weight,<u,v>>
        PriorityQueue<Pair<Integer,Pair<Integer,Integer>>> pq = new PriorityQueue<Pair<Integer,Pair<Integer,Integer>>>();
        for(int i =0 ;i < NUM_OF_VERTICES; i++) {
            for(Pair<Integer,Integer> pair: this.getAdjacentVerticesWithWeights(i)) {
                int v = pair.getFirst();
                int w = pair.getSecond();
                pq.add(new Pair<>(w, new Pair<>(i, v)));
            }
        }

        // create a disjoint set
        DisjointSet ds  = new DisjointSet();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            ds.makeSet(i);
        }

        // for each edge, if the two vertices are not in the same set, add the edge to the MST

        while(!pq.isEmpty()) {
            Pair<Integer,Pair<Integer,Integer>> pair = pq.remove();
            int u = pair.getSecond().getFirst();
            int v = pair.getSecond().getSecond();
            if(!ds.isSameSet(u, v)) {
                System.out.println(u + "->" + v);
                ds.union(u, v);
            }
        }
        return;
    } 

    public void findBrighesUtilDFS(int i,int parent, ArrayList<Boolean> isVisited,ArrayList<Integer> tin, ArrayList<Integer>low, int timer,ArrayList<Pair<Integer,Integer>> brighes) {
        isVisited.set(i,true);
        tin.set(i,timer+1);
        low.set(i,timer+1);

        ArrayList<Pair<Integer,Integer>> adjacentVertices = this.getAdjacentVerticesWithWeights(i);

        for(Pair<Integer,Integer> pair: adjacentVertices) {
            int v = pair.getFirst();
            if (v == parent) continue;
            if(!isVisited.get(v)) {
                findBrighesUtilDFS(v,i, isVisited, tin, low, timer,brighes);
                low.set(i, Math.min(low.get(i), low.get(v)));
                if (low.get(v) > tin.get(i)) {
                    brighes.add(new Pair<>(i,v));
                }
            } else {
                low.set(i, Math.min(low.get(i), tin.get(v)));
            }
        }
        return;
    }
    // find brighes in graph
    public ArrayList<Pair<Integer, Integer>> findBrighes() {
        //result array
        ArrayList<Pair<Integer, Integer>> brighes = new ArrayList<>();
        ArrayList<Boolean> isVisited = new ArrayList<>();
        ArrayList<Integer> tin = new ArrayList<>();
        ArrayList<Integer> low = new ArrayList<>();
        int timer = 0 ;
        for(int i =0 ;i < NUM_OF_VERTICES; i++)
            isVisited.add(false);
        
        for(int i = 0 ; i < NUM_OF_VERTICES; i++) {
            if(!isVisited.get(i)) {
                findBrighesUtilDFS(i,-1, isVisited, tin, low, timer,brighes);
            }
        }
        return brighes;
    }

    // articulation util
    public void findArticulationPointsUtil(int node, int parent,int timer, ArrayList<Boolean> isVisited, ArrayList<Integer> tin , ArrayList<Integer> low  , ArrayList<Integer> ap) {
        timer++;
        isVisited.set(node,true);
        tin.set(node,timer);
        low.set(node,timer);

        ArrayList<Integer> adjacentVertices = this.getAdjacentVertices(node);
        int childCount  = 0 ;

        for(int v : adjacentVertices) {
            if (v == parent) continue; 

            if (!isVisited.get(node)) {
                findArticulationPointsUtil(node, parent, timer, isVisited, tin, low, ap);
                low.set(node,Math.min(low.get(node),low.get(v)));
                // parent!=-1 means removing root won't results in a new articulation point
                if (low.get(v) >= tin.get(node) && parent != -1 ) {
                    ap.set(node,1);
                }
                childCount++;
            }
            else {
                low.set(node,Math.min(low.get(node),tin.get(v)));
            }
        }
        // handle root node that has no parent and has individual children
        // removing this results into more than one components of the graph
        if (parent == -1 && childCount > 1) {
            ap.set(node,1);
        }
        return;
    }


    // find articulations points in graph
    public ArrayList<Integer> findArticulationPoints() {
        ArrayList<Boolean> isVisited  = new ArrayList<>();
        ArrayList<Integer> tin = new ArrayList<>();
        ArrayList<Integer> low = new ArrayList<>();
        ArrayList<Integer> ap = new ArrayList<>();
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            isVisited.add(false);
            ap.add(0);
        }
        int timer = 0 ;
        for(int i = 0; i < NUM_OF_VERTICES; i++) {
            if(!isVisited.get(i)) {
                findArticulationPointsUtil(i,-1,timer, isVisited, tin, low, ap);
            }
        }
        return ap;
    }

    public void printSCCUtilDFS(int node, ArrayList<Boolean> isVisited) {
        isVisited.set(node,true);
        System.out.print(node + " ");
        for(int v : this.getAdjacentVertices(node)) {
            if (!isVisited.get(v)) {
                printSCCUtilDFS(v, isVisited);
            }
        }
    }
    // print SCC using kosaraju's algorithm
    // TC: O(V+E)
    // SC: O(V+E)
    public void printSCC() {
        /**
         * steps:
         * 1. topsort
         * 2. transpose graph
         * 3. dfs on the transposed graph
         */

        // 1. topsort
        ArrayList<Integer> topsort = this.topologicalSortUsingDFS();
        // 2. transpose graph
        Graph g = this.transposeGraph(false);
        // 3. dfs on the transposed graph
        ArrayList<Boolean> isVisited = new ArrayList<>();
        for (int i = 0; i < NUM_OF_VERTICES; i++) {
            isVisited.add(false);
        }
        for (int i =  0; i < topsort.size(); i++) {
            int node = topsort.get(i);
            if (!isVisited.get(node)) {
                g.printSCCUtilDFS(node, isVisited);
            }
        }
    }

    // implement trajan's algorithm
    
}