import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from matplotlib import pyplot as plt

def pathFinder(source,target,network):
	paths = []
	hasPath = nx.has_path(network,source,target)
	print('Is there a path from %s to %s ?\n'%(source,target))
	print(hasPath)
	if hasPath==True:
		paths = [p for p in nx.all_shortest_paths(net, source, target)]
	return paths
		


print('Welcome to target-source path founder script')
mode="wrong"
print('Enter mode customization: xlsx File of multiple source inputs (F) or one source (S)\n')
while ((mode!='F') and (mode!='S')):
	mode=input("Enter F or S: ")
	if ((mode!='F') and (mode!='S')):
		print('Wrong input value,enter again')

# Import omnipath preprocessed
df_net = pd.read_csv('preprocessed_data/FilteredOmnipath.tsv', sep='\t',index_col=0)
net = nx.from_pandas_edgelist(df_net,
                              source = 'source',
                              target = 'target',
                              edge_attr = True,
                              create_using=nx.MultiDiGraph())

outputFile = input("Enter the path/name of the output CSV file of the results (e.g. ../results/file.csv):")
if mode=='S':
	source = input("Enter gene name of source: ")
	target = input("Enter gene name of target: ")
	paths = pathFinder(source,target,net)

	for i,p in enumerate(paths):
		Path('../results/geneCandidatesSubnets/'+source).mkdir(parents=True, exist_ok=True)
		subnet = net.subgraph(p)
		# Plot the graph
		fig = plt.figure(num=None, figsize=(6, 6), dpi=80)
		ax = fig.add_subplot(111)
		pos = nx.spring_layout(subnet)
		nx.draw_networkx_nodes(subnet, pos, node_size=500, ax=ax)
		nx.draw_networkx_labels(subnet, pos,font_size=8, ax=ax)
		for edge in subnet.edges(data=True):
			sign = edge[2]['interaction']
			if sign < 0:
				sign = '-['
			else:
				sign = '-|>'
			nx.draw_networkx_edges(subnet, pos, edgelist=[(edge[0], edge[1])],
								   node_size=500,
								   arrowstyle=sign,
								   ax=ax)
		# Turn off unnecessary x and y axis labels
		ax.set_axis_off()
		# Save figure with a tight bbox to ensure it isn't cut off
		fig.savefig('../results/geneCandidatesSubnets/'+source+'/subnet_%s.png'%i, bbox_inches='tight', dpi=600)
	# Dataframe to save the result in a csv
	df = pd.DataFrame({'source': [source],
					   'paths': paths})

	df.to_csv(outputFile)
else:
	sources_file = input("Enter the path for the xlsx file of sources: e.g. ../results/my_file.xlsx")
	target = input("Enter gene name of target: ")
	sources = pd.read_excel(sources_file,index_col=0)
	sources = list(sources['genes'])
	paths = []
	for source in sources:
		path = pathFinder(source, target, net)
		paths.append(path)
		for i,p in enumerate(path):
			Path('../results/geneCandidatesSubnets/' + source).mkdir(parents=True, exist_ok=True)
			subnet = net.subgraph(p)
			# Plot the graph
			fig = plt.figure(num=None, figsize=(6, 6), dpi=80)
			ax = fig.add_subplot(111)
			pos = nx.spring_layout(subnet)
			nx.draw_networkx_nodes(subnet, pos, node_size=500, ax=ax)
			nx.draw_networkx_labels(subnet, pos,font_size=8, ax=ax)
			for edge in subnet.edges(data=True):
				sign = edge[2]['interaction']
				if sign < 0:
					sign = '-['
				else:
					sign = '-|>'
				nx.draw_networkx_edges(subnet, pos, edgelist=[(edge[0], edge[1])],
									   node_size=500,
									   arrowstyle=sign,
									   ax=ax)
			# Turn off unnecessary x and y axis labels
			ax.set_axis_off()
			# Save figure with a tight bbox to ensure it isn't cut off
			fig.savefig('../results/geneCandidatesSubnets/' + source + '/subnet_%s.png'%i, bbox_inches='tight', dpi=600)
	# Dataframe to save the result in a csv
	df = pd.DataFrame({'source':sources,
					   'paths':paths})
	df.to_csv(outputFile)
