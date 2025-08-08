#!/usr/bin/env python3
import sys
import xmltodict
from networkx import nx
import matplotlib.pyplot as plt
from LinkActiveV import *
import collections
import numpy as np
from matplotlib import pyplot

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
	
	degree_sequence=sorted(nx.degree(G).values(),reverse=True) # degree sequence
	#print "Degree sequence", degree_sequence
	degreeCount=collections.Counter(degree_sequence)
	deg, cnt = zip(*degreeCount.items())
	fig, ax1 = plt.subplots()
	plt.bar(deg, cnt, width=0.80, color='b')
	plt.title("Degree Histogram of Full Network")
	plt.ylabel("Count")
	plt.xlabel("Degree")
	ax1.set_xticks([d+0.4 for d in deg])
	ax1.set_xticklabels(deg)
	# draw graph in inset
	#plt.axes([0.4, 0.4, 0.5, 0.5])
	Gc=sorted(nx.connected_component_subgraphs(G), key = len, reverse=True)[0]
	pos=nx.spring_layout(G)
	#plt.axis('off')
	#nx.draw_networkx_nodes(Gc, pos, node_size=20)
	#nx.draw_networkx_edges(Gc, pos, alpha=0.4)
	plt.savefig("degree_histogram_1.png")
	plt.show()
	
	Ga = nx.Graph()
	Ga.add_edges_from(ActiveEdges)
	pos_a = nx.spring_layout(Ga)
	valuesA = [val_map.get(node, 0.65) for node in Ga.nodes()]
	
	degree_sequence=sorted(nx.degree(Ga).values(),reverse=True) # degree sequence
	degreeCount=collections.Counter(degree_sequence)
	deg, cnt = zip(*degreeCount.items())

	fig, ax2 = plt.subplots()
	plt.bar(deg, cnt, width=0.80, color='b')
	plt.title("Degree Histogram of Active Network")
	plt.ylabel("Count")
	plt.xlabel("Degree")
	ax2.set_xticks([d+0.4 for d in deg])
	ax2.set_xticklabels(deg)

	Gcc=sorted(nx.connected_component_subgraphs(Ga), key = len, reverse=True)[0]
	#plt.axis('off')
	#nx.draw_networkx_nodes(Ga, pos_a, node_size=250)
	#nx.draw_networkx_edges(Ga, pos_a, alpha=0.9)
	plt.savefig("degree_histogram_2.png")
	plt.show()

	#print(nx.clustering(G))
	#print(nx.clustering(Ga))
	
	
def main(argv):

	net_plot(argv[1],argv[2])
	
if __name__ == '__main__':
    sys.exit(main(sys.argv))

