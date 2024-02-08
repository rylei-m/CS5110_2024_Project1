"""
Written by Vicki Allan using Bellman Ford pseudocode
For me, it is easier to write the min cost max flow using adjacency matrices, so that is what this code does
Feel free to do it your way.
I was learning Python as I wrote, so I'm sure you'll shake your head at some of the things I did.  Shaking
your head is good exercise!
"""

##IMPORTS ARE USED TO SHOW THE GRAPH AT THE END
import networkx as nx
import matplotlib.pyplot as plt

class Graph:

    # We should have two edges: one from a->b indicatge preference for men and women.
    # If we don't have two edges, one partner found the other unacceptable, so we can ignore that combination
    def combine(self, edge1, edge2):
        if (edge1[0] == edge2[0] and edge1[1] == edge2[1]):
            return (self.vertices.index(edge1[0]), self.vertices.index(edge1[1]), edge1[2] + edge2[2], 1)
        else:
            return None

    # for our matching, we consider the cost of an edge to be the sum of costs for each partner.
    # In this method, we combine individual edges creating a new set of edges.
    def combine_edges(self, edges):
        new_edges = []
        i = 0

        while i < len(edges) -1:
            e = self.combine(edges[i], edges[i + 1])
            if e is not None:
                new_edges.append(e)
                i += 2
            else:
                i += 1  # skip the edge without a match
        return new_edges

    def create_graph(self, file_tuple):
        self.vertices.append("Source")  # create a list of all vertices
        # f = file(filename)
        with open(file_tuple[0]) as f:
            for line in f:
                pieces = line.split(':')
                name = pieces[0].strip()
                self.vertices.append(name)

                if name:
                    priorities = pieces[1].strip().split(',')
                    for i in range(len(priorities)):
                        priorities[i] = priorities[i].strip()
                        # create an edge a->b with cost and flow
                        self.edges.append((name, priorities[i], i + 1, 1))
            f.close()
        men_ct = len(self.vertices) - 1
        # repeat for women's preferences
        with open(file_tuple[1]) as f:
            for line in f:
                pieces = line.split(':')
                name = pieces[0].strip()
                self.vertices.append(name)

                if name:
                    priorities = pieces[1].strip().split(',')
                    for i in range(len(priorities)):
                        priorities[i] = priorities[i].strip()
                        # create an edge a->b with cost and flow
                        self.edges.append((priorities[i], name, i + 1, 1))
            f.close()
        self.vertices.append("Sink")
        # Sort edges to get those with same to/from nodes together
        self.edges.sort()
        self.edges = self.combine_edges(self.edges)
        sink = len(self.vertices) - 1

        # set up edges as max flow problem
        for men in range(1, men_ct + 1):
            self.edges.append((0, men, 0, 1))
        for women in range(men_ct + 1, sink):
            self.edges.append((women, sink, 0, 1))
        self.make_adjacency()

    # from the list of edges, create an adjacency matrix, residual matrix, and cost_matrix
    def make_adjacency(self):
        self.vertex_ct = len(self.vertices)
        self.adjM = []
        while (len(self.adjM) < self.vertex_ct):
            temp = [0 for i in range(self.vertex_ct)]
            self.adjM.append(temp)
        self.cost_array = [list(row) for row in self.adjM]  # careful to get a deep copy

        for edge in self.edges:
            i = int(edge[0])
            j = int(edge[1])

            if i >= self.vertex_ct or j >= self.vertex_ct or i < 0 or j < 0:
                print(f"Not a Proper Input in Edge {i},{j}")
            else:
                self.adjM[i][j] = edge[3]
                self.cost_array[i][j] = edge[2]
                self.cost_array[j][i] = -edge[2]
            self.residual = [list(row) for row in self.adjM]  # careful to get a deep copy

    # print 2 d array a with label
    @staticmethod
    def print2d_array(label, a):
        print(label)
        for i in range(len(a)):
            print("%3d:" % (i), end=' ')
            for j in range(len(a[i])):
                print("%3d" % (a[i][j]), end=' ')
            print()

    def do_flow(self):
        print("Vertices are: ")
        print(self.vertices)
        print("Edges are: ")
        print(self.edges)
        self.print2d_array("adjacency", self.adjM)
        self.print2d_array("residual", self.residual)
        self.print2d_array("cost", self.cost_array)
        self.BellmanFord(0, len(self.vertices) - 1)

    # utility function used to print the matrix dist with label
    def print_array(self, label, dist):
        print(label)
        for i in range(self.vertex_ct):
            print("{0}\t\t{1}".format(i, dist[i]))

    def visualize_graph(self): ##VISUALIZE GRAPH IS FROM OPENAI SO THAT THE GRAPH CAN BE SHOWN
        G = nx.DiGraph()

        # Add nodes
        for node in self.vertices:
            G.add_node(node)

        # Add edges
        for edge in self.edges:
            u, v, weight, _ = edge
            G.add_edge(self.vertices[u], self.vertices[v], weight=weight)

        pos = nx.spring_layout(G)  # positions for all nodes
        
        # nodes
        nx.draw_networkx_nodes(G, pos, node_size=700)

        # edges
        nx.draw_networkx_edges(G, pos, edgelist=G.edges(), width=2)
        nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): d['weight'] for u, v, d in G.edges(data=True)})

        # labels
        nx.draw_networkx_labels(G, pos, font_size=20, font_family='sans-serif')

        plt.axis('off')
        plt.show()

    # The main function that finds shortest distances from src to
    # all other vertices using Bellman-Ford algorithm. The function
    # also detects negative weight cycle
    # We are interested in the path from the src to the sink.
    # If we never make it to the sink, there is no flow
    # return true if there is flow from src to sink.
    def BellmanFord(self, src, sink):

        # Step 1: Initialize distances from src to all other vertices
        # as INFINITE
        dist = [9999 for i in range(self.vertex_ct)]  # dist/cost to each node
        pred = [-1 for i in range(self.vertex_ct)]  # predecessor of each node
        dist[src] = 0

        # Step 2: Relax all edges |V| - 1 times. A simple shortest
        # path from src to any other vertex can have at-most |V| - 1
        # edges    
        for _ in range(self.vertex_ct - 1):
            for u in range(self.vertex_ct):
                for v in range(self.vertex_ct):
                    if self.residual[u][v] > 0 and dist[u] != 9999 and dist[u] + self.cost_array[u][v] < dist[v]:
                        dist[v] = dist[u] + self.cost_array[u][v]
                        pred[v] = u
         
        path = []
        if pred[sink] != -1:  # If there's a path to the sink
            current = sink
            while current != src:
                path.insert(0, current)
                current = pred[current]
            path.insert(0, src)
        return path
    
        self.print_array("Predecessor", pred)
        self.print_array("Cost", dist)
        return pred[sink] >= 0

    def min_cost_max_flow(self):
        self.flow = [[0 for _ in range(self.vertex_ct)] for _ in range(self.vertex_ct)]  # Initialize flow matrix
        total_flow = 0
        total_cost = 0
        while True:
            # Find path with Bellman-Ford
            path = self.BellmanFord(0, len(self.vertices) -1)
            if not path or len(path) ==1:
                break

            # Find the minimum residual capacity of the path found
            flow = self.min_capacity_of_path(path)
            path_cost = self.calculate_path_cost(path)
            self.augment_flow_along_path(path, flow)

            # Augment flow and update residual capacities and costs
            total_flow += flow
            total_cost += flow * path_cost

        return total_flow, total_cost
    
        # create an adjacency matrix from men preferencese and women preferences
    def __init__(self, fileTuple):
        self.vertices = []
        self.adjM = []
        self.vertex_ct = 0
        self.edges = []
        self.residual = []
        self.cost_array = []
        self.flow = []  # Initialize flow matrix
        self.create_graph(fileTuple)
    
    def min_capacity_of_path(self, path):
    # Initialize minimum capacity to a large number
        min_capacity = float('inf')

        # Iterate through the path
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]

            # Find the residual capacity of the edge (u, v)
            residual_capacity = self.residual[u][v]

            # Update min_capacity if this edge has a lower residual capacity
            if residual_capacity < min_capacity:
                min_capacity = residual_capacity

        return min_capacity
    
    def augment_flow_along_path(self, path, flow):
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]

            # Augment the flow along the edge (u, v)
            # Assuming there's a structure to hold the current flow, update it
            self.flow[u][v] += flow
            self.flow[v][u] -= flow  # Update the reverse flow for undirected graph

            # Update the residual capacities
            self.residual[u][v] -= flow
            self.residual[v][u] += flow

    def calculate_path_cost(self, path):
        total_cost = 0
        for i in range(len(path) - 1):
            u = path[i]
            v = path[i + 1]

            # Add the cost of the edge (u, v) to the total cost
            # Assuming self.cost_array stores the cost of each edge
            total_cost += self.cost_array[u][v]

        return total_cost
    
    def FordFulkerson(self, source, sink):
        total_flow = 0
        while True:
            path = self.BellmanFord(source, sink)
            if not path or len(path) == 1:  # No augmenting path found
                break

            # Find the bottleneck capacity
            avail_flow = float('inf')
            for v in range(1, len(path)):
                u = path[v - 1]
                v = path[v]
                avail_flow = min(avail_flow, self.residual[u][v])

            # Augment flow and update residual capacities
            for v in range(1, len(path)):
                u = path[v - 1]
                v = path[v]
                self.residual[u][v] -= avail_flow
                self.residual[v][u] += avail_flow

            total_flow += avail_flow

        return total_flow

##files = [("men0.txt", "women0.txt", True), ("men.txt", "women.txt", True), ("men2.txt", "women2.txt", True),
##         ("men3.txt", "women3.txt", False), ("men4.txt", "women4.txt", False)]
files = [("Applicants2.txt", "Employers2.txt", True)]
for fileTuple in files:
    g=Graph(fileTuple)
    g.do_flow()


def stable_matching_employers_propose(employers_prefs, applicants_prefs):
    # Initialize all employers as unmatched
    matches = {employer: None for employer in employers_prefs}
    unmatched_employers = set(employers_prefs.keys())

    while unmatched_employers:
        for employer in list(unmatched_employers):
            if not employers_prefs[employer]:
                unmatched_employers.remove(employer)
                continue  # Skip if employer has no more preferred applicants

            applicant = employers_prefs[employer].pop(0)
            current_match = next((e for e, a in matches.items() if a == applicant), None)

            if not current_match:
                matches[employer] = applicant
                unmatched_employers.remove(employer)
            else:
                # Ensure that the current employer and current match are in applicant's preferences
                if applicant in applicants_prefs and employer in applicants_prefs[applicant] and current_match in applicants_prefs[applicant]:
                    if applicants_prefs[applicant].index(employer) < applicants_prefs[applicant].index(current_match):
                        matches[employer] = applicant
                        unmatched_employers.remove(employer)
                        unmatched_employers.add(current_match)
                        matches[current_match] = None
                else:
                    unmatched_employers.remove(employer)
    return matches

def stable_matching_applicants_propose(applicants_prefs, employers_prefs):
    matches = {applicant: None for applicant in applicants_prefs}
    unmatched_applicants = set(applicants_prefs.keys())

    while unmatched_applicants:
        for applicant in list(unmatched_applicants):
            if not applicants_prefs[applicant]:
                unmatched_applicants.remove(applicant)
                continue  # Skip if applicant has no more preferred employers

            employer = applicants_prefs[applicant].pop(0)
            current_match = next((a for a, e in matches.items() if e == employer), None)

            if not current_match:
                matches[applicant] = employer
                unmatched_applicants.remove(applicant)
            else:
                if employer in employers_prefs and applicant in employers_prefs[employer] and current_match in employers_prefs[employer]:
                    if employers_prefs[employer].index(applicant) < employers_prefs[employer].index(current_match):
                        matches[applicant] = employer
                        unmatched_applicants.remove(applicant)
                        unmatched_applicants.add(current_match)
                        matches[current_match] = None
                else:
                    unmatched_applicants.remove(applicant)
    return matches

def min_cost_max_flow_matching(fileTuple):
        g = Graph(fileTuple)
        total_flow, total_cost = g.min_cost_max_flow()
        print(f"Results for {fileTuple}: Total Flow = {total_flow}, Total Cost = {total_cost}")
    
    ## PARSING
def parse_preferences(file_path):
    preferences = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(':')
            if len(parts) == 2:
                key = parts[0].strip()
                values = parts[1].strip().split(',')
                preferences[key] = [value.strip() for value in values]
    return preferences

# Replace 'your_applicants_file_path' and 'your_employers_file_path' with the actual paths to your files
applicants_file_path = 'Applicants2.txt'
employers_file_path = 'Employers2.txt'

applicants_preferences = parse_preferences(applicants_file_path)
employers_preferences = parse_preferences(employers_file_path)
stable_matches = stable_matching_employers_propose(employers_preferences, applicants_preferences)
stable_matches_applicants_proposing = stable_matching_applicants_propose(applicants_preferences, employers_preferences)
print("Stable Matching with Applicants Proposing:")
print(stable_matches_applicants_proposing)
print("Stable Matching with Employers Proposing:", stable_matches)
print(stable_matches)
print("applicants preferences:", applicants_preferences)
print("employers preferences:", employers_preferences)
min_cost_max_flow_matching((applicants_file_path, employers_file_path))

graph = Graph((applicants_file_path, employers_file_path))
graph.visualize_graph()
