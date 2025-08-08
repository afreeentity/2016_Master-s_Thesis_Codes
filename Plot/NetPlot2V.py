#!/usr/bin/env python3
import sys
import xmltodict
from networkx import nx
import matplotlib.pyplot as plt
from LinkActiveV import *
import numpy as np

def color_node(a,b,name):
	if a != '0':
		x = a
	else:
		x = b	
	
	if name == 'Vt':
		y = 0.05
	else:
		y = 0.9

	return (x,y)
		

def net_plot(ff,fs):

	thr = float(input("Please give me the threshold of active nodes: "))

	ActiveEdges = []
	with open(fs,'r') as fstart:
		for line in fstart:
			kl = line.split()
			if link_active(float(kl[3]),thr): ActiveEdges.append((kl[1],kl[2]))


	with open(ff, 'r') as Input:
		doc = xmltodict.parse(Input.read())
	
	Elem = doc['mydocument']['Elements']	
	elist=[e for e in Elem]
	print(elist)
	val_map = {}
	node_list = []



	for name in elist:
		if isinstance(Elem[name],list):
			for i in Elem[name]:
				if i['G'][0]['#text']=='0' or i['G'][1]['#text']=='0' :
					k,h = color_node( i['G'][0]['#text'], i['G'][1]['#text'],name)					
					val_map[k] = h					
				else:
					node_list.append( (i['G'][0]['#text'], i['G'][1]['#text']) )
		else:
			i=Elem[name]
			if i['G'][0]['#text']=='0' or i['G'][1]['#text']=='0' :
				k,h = color_node( i['G'][0]['#text'], i['G'][1]['#text'],name)					
				val_map[k] = h					
			else:
				node_list.append( (i['G'][0]['#text'], i['G'][1]['#text']) )

	G = nx.Graph()
	G.add_edges_from(node_list)

	pos = nx.spring_layout(G)
	values = [val_map.get(node, 0.65) for node in G.nodes()]
	
	degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
	#print "Degree sequence", degree_sequence
	dmax=max(degree_sequence)

	plt.loglog(degree_sequence,'b-',marker='o')
	plt.title("Degree Rank Plot of Full Network")
	plt.ylabel("degree")
	plt.xlabel("rank")

	Gc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
	pos=nx.spring_layout(Gc)
	#plt.axis('off')
	#nx.draw_networkx_nodes(Gc,pos,node_size=250)
	#nx.draw_networkx_edges(Gc,pos,alpha=0.9)

	plt.savefig("degree_1.png")
	plt.show()
	
	
	Ga = nx.Graph()
	Ga.add_edges_from(ActiveEdges)
	pos_a = nx.spring_layout(Ga)
	valuesA = [val_map.get(node, 0.65) for node in Ga.nodes()]
	
	degree_sequence=sorted(nx.degree(Ga).values(),reverse=True) # degree sequence
	#print "Degree sequence", degree_sequence
	dmax=max(degree_sequence)

	plt.loglog(degree_sequence,'b-',marker='o')
	plt.title("Degree Rank Plot of Active Network")
	plt.ylabel("degree")
	plt.xlabel("rank")

	Gcc=sorted(nx.connected_component_subgraphs(Ga), key = len, reverse=True)[0]
	pos_a=nx.spring_layout(Gcc)
	#plt.axis('off')
	#nx.draw_networkx_nodes(Gcc,pos_a,node_size=250)
	#nx.draw_networkx_edges(Gcc,pos_a,alpha=0.9)

	plt.savefig("degree_2.png")
	plt.show()
	
	

	#print(nx.clustering(G))
	#print(nx.clustering(Ga))
def main(argv):

	net_plot(argv[1],argv[2])
	
if __name__ == '__main__':
    sys.exit(main(sys.argv))

