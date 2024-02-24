"""
    This module is your primary workspace. Add whatever helper functions, classes, data structures, imports... etc here.

    We expect most results will utilize more than just dumping code into the plan_paths()
        function, that just serves as a meaningful entry point.

    In order for the rest of the scoring to work, you need to make sure you have correctly
        populated the Destination.path for each result you produce.
"""
import typing
from queue import PriorityQueue
import yaml
import numpy as np
from typing import Dict
import argparse
from map_info import Coordinate, Destination, MapInfo
import heapq

class PathPlanner:
    def __init__(self, map_info: MapInfo, destinations: typing.List["Destination"]):
        self.map_info: MapInfo = map_info
        self.destinations: typing.List["Destination"] = destinations


    def plan_paths(self):
        """
        This is the function you should re-write. It is expected to mutate the list of
        destinations by calling each Destination's set_path() with the resulting
        path as an argument.

        The default construction shows this format, and should produce 10 invalid paths.
        """

        #setup
        graph = create_graph(self.map_info)
        riskPaths = graph.dijkstras(graph.nodes[self.map_info.start_coord])
        SearchNode.lengthimp = True
        lengthPaths = graph.dijkstras(graph.nodes[self.map_info.start_coord])
        SearchNode.lengthimp = False
        
        for site in self.destinations:
            
            path_array = []

            if graph.nodes[site.coord] not in riskPaths: # node was not reached by risk first search
                if graph.nodes[site.coord] not in lengthPaths: # node was not reached at all
                    continue
                search = lengthPaths[graph.nodes[site.coord]]
            else: # node was reached by risk first search
                search = riskPaths[graph.nodes[site.coord]]
            
            while search.predecessor: # form path into an array

                search.node.coord
                path_array.append(search.node.coord)
                search = search.predecessor
            
            path_array.append(search.node.coord)
            path_array.reverse()

            site.set_path(path_array)


# This method creates a graph representing the field, creates a node for each coordinate and an edge between each adjacent node

def create_graph(map_info):
    graph = Graph(map_info)
    for i in range(0,59):
        for j in range(0,40):
             
            graph.add_node(Coordinate(i, j))
            if i != 60:
                graph.add_edge(Coordinate(i,j), Coordinate(i+1,j), 1)
                if j != 0:
                    graph.add_edge(Coordinate(i,j-1), Coordinate(i-1,j), 1.414)
            if j != 0:
                graph.add_edge(Coordinate(i,j-1), Coordinate(i,j), 1)
                graph.add_edge(Coordinate(i,j-1), Coordinate(i+1,j), 1.414)
        
    return graph
    

# Copied method to find the risk of a path to be used for finding edge weights

def find_path_risk(path_coords, map_info) -> int:
    total_risk = sum(
        map_info.risk_zones[int(coord[0])][int(coord[1])] == MapInfo.HIGH_RISK_VALUE
        for coord in path_coords
    )
    return total_risk


# Copied method to find whether or not an edge is valid to be used for determining whether or not to put an edge in the graph

def is_keepout_valid(map_info: MapInfo, path: typing.List["tuple"]) -> bool:

    return all(
        [map_info.risk_zones[int(coord[0])][int(coord[1])] < MapInfo.KEEP_OUT_VALUE for coord in path]
    )

# Copied method to find whether or not an edge is valid to be used for determining whether or not to put an edge in the graph

def is_keepin_valid(map_info: MapInfo, path: typing.List["tuple"]) -> bool:
    return all(
        [
            all([coord[0] >= 0 for coord in path]),
            all([coord[0] <= map_info.risk_zones.shape[0] for coord in path]),
            all([coord[1] >= 0 for coord in path]),
            all([coord[1] <= map_info.risk_zones.shape[1] for coord in path]),
        ]
    )

# This class represents a search node object for dijkstras algorithms. The class also contains a static boolean lengthimp which
# Tells whether or not length is more important than risk

class SearchNode:

    lengthimp = False # This represents whether we are more concerned about the weight comparison or the length comparison

    def __init__(self, node, weight, length, predecessor):
        self.node = node
        self.weight = weight
        self.length = length
        self.predecessor = predecessor
    

    def __eq__(self, other):
        return (self.weight == other.weight and self.length == other.length)

    def __lt__(self, other):
        if (SearchNode.lengthimp):
            if (self.length == other.length):
                return self.weight<other.weight
            return self.length < other.length
        if (self.weight == other.weight):
            return self.length < other.length
        return self.weight < other.weight
    
    def __gt__(self, other):
        if (SearchNode.lengthimp):
            if (self.length == other.length):
                return self.weight>other.weight
            return self.length > other.length
        if (self.weight == other.weight):
            return self.length > other.length
        return self.weight > other.weight


# This class represents an edge for our graph. Each edge has a beginning, an end, a length representing the length of the edge
# and weight,representing the risk of the edge

class Edge:

    def __init__(self, beginning, end, map_info, length):  
        self.weight = find_path_risk([beginning.coord, end.coord], map_info)
        self.beginning = beginning
        self.end = end
        self.length = length

# This class represents a node for our graph. Each node has coordinates that it corresponds to and an adjacency list of Edges

class GraphNode:

    def __init__(self, coord):
        self.coord = coord
        self.adjacency_list = []

# This class represents the field that we are asked to make paths on. Relatively standard Graph ADT implementation

class Graph:

    def __init__(self, map_info):
        self.nodes = {}
        self.map_info = map_info
        

    # This is an implementation of Dijkstras shortest path implementation to find the shortest path based on risk
    # or length depending on wether or not lengthimp in the edge class is True to each point in our graph

    def dijkstras(self, startNode) -> dict:
        pq = []
        visitedNodes = {}
        heapq.heappush(pq, SearchNode(startNode, 0, 0, None))
        if (SearchNode.lengthimp): visitedNodes = {}
        while len(pq) != 0:

            curr = heapq.heappop(pq)
            if curr.node not in visitedNodes:
                visitedNodes[curr.node] = curr
                for edge in curr.node.adjacency_list:

                    if edge.end not in visitedNodes:

                        if (curr.length+edge.length <= self.map_info.maximum_range or SearchNode.lengthimp):
                            heapq.heappush(pq, SearchNode(edge.end, curr.weight+edge.weight, curr.length+edge.length, curr))
        return visitedNodes

    def add_node(self, coord:Coordinate):
        if coord not in self.nodes:
            self.nodes[coord] = GraphNode(coord)

    def add_edge(self, coord1, coord2, length):
        # check if edge will be valid invalid edges simply won't be placed
        if (not is_keepout_valid(self.map_info, [coord1, coord2])): return
        if (not is_keepin_valid(self.map_info, [coord1, coord2])): return
        #edge is valid, add to graph
        if coord1 in self.nodes:
            if coord2 not in self.nodes:
                self.add_node(coord2)
            self.nodes[coord1].adjacency_list.append(Edge(self.nodes[coord1], self.nodes[coord2], self.map_info, length))
            self.nodes[coord2].adjacency_list.append(Edge(self.nodes[coord2], self.nodes[coord1], self.map_info, length))

        else:
            raise ValueError("Start Node not in graph")

