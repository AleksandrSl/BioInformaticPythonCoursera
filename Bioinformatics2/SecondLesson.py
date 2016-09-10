import random, Bioinf2.FirstLesson as l2
from copy import deepcopy
from itertools import product
adjacency_list = """""".splitlines()

def make_graph(adjacency_list):
    graph_dict = {}
    for line in adjacency_list:
        nodes = list(map(str.strip, [line.split('->')[0]] + line.split('->')[1].split(',')))
        for node in nodes:
            if node not in graph_dict:
                graph_dict[node] = {"in":[],"out":[]}
        for node in nodes[1:]:
            graph_dict[nodes[0]]["out"].append(node)
            graph_dict[node]["in"].append(nodes[0])
    return graph_dict
graph = make_graph(adjacency_list)


def eulerian_cycle(graph):
    eul_cycle = []
    node = random.choice(list(graph.keys()))
    eul_cycle.append(node)
    insert_point = len(eul_cycle)
    while graph != {}:
        while node in graph:
            next_node = graph[node]["out"].pop()
            if not graph[node]["out"]:          # if graph[node]["out"] == 0
                graph.pop(node)
            node = next_node
            eul_cycle.insert(insert_point, node)
            insert_point += 1
        for i in range(len(eul_cycle)):
            if eul_cycle[i] in graph:
                node = eul_cycle[i]
                insert_point = i + 1
                break
    return eul_cycle



def eulerian_path(graph):
    eul_cycle = []
    not_visited = deepcopy(graph)
    print(graph)

    next_node = None
    for item in graph.items():
        if not item[1]["out"]:
            break_point = item[0]
            not_visited.pop(item[0])
    print(break_point)
    print(not_visited)
    node = random.choice(list(not_visited.keys()))
    eul_cycle.append(node)
    insert_point = len(eul_cycle)
    while not_visited != {}:
        if node == next_node:
            node = random.choice(list(not_visited.keys()))
            insert_point = len(eul_cycle)
        while node in not_visited:
            next_node = not_visited[node]["out"].pop()
            if not not_visited[node]["out"]:  # if graph[node]["out"] == 0
                not_visited.pop(node)
            node = next_node
            eul_cycle.insert(insert_point, node)
            insert_point += 1
        for i in range(len(eul_cycle)):
            if eul_cycle[i] in not_visited:
                node = eul_cycle[i]
                insert_point = i + 1
                break
        print(eul_cycle)
    break_point = eul_cycle.index(break_point) + 1

    return eul_cycle[break_point:] + eul_cycle[:break_point ]




def eulerian_path2(graph):
    eul_cycle = []
    stack = []
    not_visited = deepcopy(graph)
    #print(graph)
    node = None
    start_node = None
    for item in graph.items():
        if (len(item[1]["out"]) - len(item[1]["in"])) == 1:
            node = item[0]
    if node is None:
        node = random.choice(list(not_visited.keys()))
    print(node)
    stack.append(node)
    while not_visited[node]["out"] or (len(stack) != 0):
        if not not_visited[node]["out"]:
            eul_cycle.append(node)
            node = stack.pop()
            #print(node)
        else:
            stack.append(node)
            # print(node, not_visited[node]["out"])
            node = not_visited[node]["out"].pop()
    eul_cycle.reverse()
        #print(stack)
    #print(eul_cycle)
    return eul_cycle

#graph = l2.de_bruijn_graph2(reads)
#print(graph)
#print('->'.join(eulerian_path2(graph)).lstrip('->'))

def reconstr_from_reads(reads):
    reads = reads.splitlines()
    graph = l2.de_bruijn_graph2(reads)
    eulerian_path = eulerian_path2(graph)
    reconstr = ""
    for kmer in eulerian_path:
        reconstr += kmer[0]
    reconstr += eulerian_path[-1][1:]
    return reconstr

print()
print(reconstr_from_reads("""AAAT
AATG
ACCC
ACGC
ATAC
ATCA
ATGC
CAAA
CACC
CATA
CATC
CCAG
CCCA
CGCT
CTCA
GCAT
GCTC
TACG
TCAC
TCAT
TGCA"""))
#print(reconstr_from_reads(read))
def hammingDist(seq1, seq2):
    assert len(seq1) == len(seq2)
    hamDist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            hamDist += 1
    return hamDist


def k_universal_circular_str(k):
    #strings = '\n'.join(neighbours_binary('0'*k, k)).lstrip()
    strings = ['\n'.join(x) for x in product('01', repeat=k)]
    graph = l2.de_bruijn_graph2(strings.splitlines())
    print(graph)
    print(eulerian_path2(graph))
    res = reconstr_from_reads(strings)
    return res[:-k + 1]

def gen_pair_reads(read, k, d):
    pair_reads = []
    for i in range(len(read) - (2 * k + d) + 1):
        pair_reads.append((read[i:i + k], read[i + k + d:i + d + 2*k]))
    pair_reads.sort()
    return pair_reads

#for pair in gen_pair_reads("TAATGCCATGGGATGTT", 3, 2):
   # print('|'.join(pair).lstrip())

def de_bruijn_graph_from_pair_reads(pair_reads):
    graph = {}
    for pair_read in pair_reads:
        reads = pair_read.split('|')
        print("reads",reads)
        suf = (reads[0][1:], reads[1][1:])
        pref = (reads[0][:-1], reads[1][:-1])
        if suf not in graph:
            graph[suf] = {'in': [], 'out': [], 'weights': []}
        if pref not in graph:
            graph[pref] = {'in': [], 'out': [], 'weights': []}
        graph[pref]['out'].append(suf)
        graph[suf]['in'].append(pref)
        graph[pref]['weights'].append(reads)
    return graph

def string_spelled_by_gapped_patterns(k, d, pair_reads):
    first_part = []
    last_part = []
    for pair_read in pair_reads:
        reads = pair_read.split("|")
        first_part.append(reads[0])
        last_part.append(reads[1])
    suf = ""
    pref = ""
    string = ""
    for kmer in last_part:
         suf += kmer[0]
    suf += last_part[-1][1:]
    for kmer in first_part:
        pref += kmer[0]
    pref += first_part[-1][1:]
    for i in range(k+d,len(pref)):
        if pref[i] != suf[i - (k + d)]:
            print("there is no string spelled by the gapped patterns")
            return
    string = pref[:k + d] + suf
    return string

def string_spelled_by_gapped_patterns2(k, d, pair_reads):
    first_part = []
    last_part = []
    for pair_read in pair_reads:
        first_part.append(pair_read[0])
        last_part.append(pair_read[1])
    suf = ""
    pref = ""
    string = ""
    for kmer in last_part:
         suf += kmer[0]
    suf += last_part[-1][1:]
    for kmer in first_part:
        pref += kmer[0]
    pref += first_part[-1][1:]
    for i in range(k+d,len(pref)):
        if pref[i] != suf[i - (k + d)]:
            print("there is no string spelled by the gapped patterns")
            return
    string = pref[:k + d] + suf
    return string



#reads = """""".splitlines()

def reconstruct_from_pair_reads(k, d, pair_reads): # Нихуя не работает
    graph = de_bruijn_graph_from_pair_reads(pair_reads)
    print(graph)
    eul_path = eulerian_path2(graph)
    print(eul_path)
    return string_spelled_by_gapped_patterns2(k, d, eul_path)


print(reconstruct_from_pair_reads(3, 1,
['ACC|ATA',
'ACT|ATT',
'ATA|TGA',
'ATT|TGA',
'CAC|GAT',
'CCG|TAC',
'CGA|ACT',
'CTG|AGC',
'CTG|TTC',
'GAA|CTT',
'GAT|CTG',
'GAT|CTG',
'TAC|GAT',
'TCT|AAG',
'TGA|GCT',
'TGA|TCT',
'TTC|GAA']))

#print(reconstruct_from_pair_reads(4,2,['GACC|GCGC',
#'ACCG|CGCC',
#'CCGA|GCCG',
#'CGAG|CCGG',
#'GAGC|CGGA']) )


def max_non_branching_paths(graph):
    paths = []
    nodes_in_cycles = set()
    print(graph)
    cycles = []
    for node in graph:
        if (len(graph[node]["in"]) != 1) or (len(graph[node]["out"]) != 1):
            #print(node)
            if graph[node]["out"]:
                #print("!",node)
                for out in graph[node]["out"]:
                    non_branch_path = [node]
                    next_n = out
                    while (len(graph[next_n]["in"]) == 1) and (len(graph[next_n]["out"]) == 1):
                        non_branch_path.append(next_n)
                        next_n = graph[next_n]["out"][0]
                    non_branch_path.append(next_n)
                    paths.append(non_branch_path)
        elif (len(graph[node]["in"]) == 1) and (len(graph[node]["out"]) == 1):
            if node not in nodes_in_cycles:
                start_n = node
                is_cycle = True
                out = graph[start_n]["out"][0]
                non_branch_path = [node]
                while out != start_n:
                    if (len(graph[out]["in"]) == 1) and (len(graph[out]["out"]) == 1):
                        non_branch_path.append(out)
                        out = graph[out]["out"][0]
                    else:
                        is_cycle = False
                        break
                if is_cycle:
                    non_branch_path.append(out)
                    m = str(min(list(map(int,non_branch_path))))
                    #print(non_branch_path)
                    #print(m)
                    m_i = non_branch_path.index(m)
                    if m_i != 0:
                        non_branch_path = non_branch_path[m_i:] + non_branch_path[1:m_i] + [m]
                    cycles.append(non_branch_path)
                    #print("!",non_branch_path)

                    nodes_in_cycles.update(non_branch_path)
    paths.extend(cycles)
    return paths

adjacency_list = """""".splitlines()

#graph = make_graph(adjacency_list)
#for i in max_non_branching_paths(graph):
#    print(" -> ".join(i).rstrip(" -> "))

def gen_contigs(reads):
    short_reads = []
    #for read in reads:
    #    short_reads.extend(l2.composition(read, max(2,len(read)//2)))
    #print(short_reads)
    graph = l2.de_bruijn_graph2(reads)
    paths = max_non_branching_paths(graph)
    for path in paths:
        reconstr = ""
        for kmer in path:
            reconstr += kmer[0]
        reconstr += path[-1][1:]
        print(reconstr)


reads = """""".splitlines()
gen_contigs(reads)
