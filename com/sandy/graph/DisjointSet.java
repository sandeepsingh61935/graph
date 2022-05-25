package com.sandy.graph;
import java.util.ArrayList;

public class DisjointSet {
    private ArrayList<Integer> rank;
    private ArrayList<Integer> parent;
    DisjointSet() {
        this.rank = new ArrayList<Integer>();
        this.parent = new ArrayList<Integer>();
    }

    public void makeSet(int x) {
        rank.add(0);
        parent.add(x);
    }

    public int getRank(int index) {
        return rank.get(index);
    }

    public void setRank(int index, int rank) {
        this.rank.set(index, rank);
        return;
    }

    public int getParent(int index) {
        return parent.get(index);
    }

    public void setParent(int index, int parent) {
        this.parent.set(index, parent);
        return;
    }


    public int find(int x) {
        if(parent.get(x)==x) {
            return x;
        }
        else {
            // path compression
            int pr = find(parent.get(x));
            setParent(x, pr);
            return pr;
        }
    }

    public void union(int x, int y) {
        x = find(x);
        y = find(y);
        int rankx = getRank(x);
        int ranky = getRank(y);
        if (rankx == ranky) {
            setParent(y, x);
            setRank(x,rankx+1);
        }
        else if(rankx > ranky) {
            setParent(y, x);
        }
        else {
            setParent(x, y);
        }
        return ;
    }

    boolean isSameSet(int x ,int y ) {
        return find(x) == find(y);
    }
}

