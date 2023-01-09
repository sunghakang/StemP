function [G,adjacency_matrix] = find_graph(edge)
% Author: Kumbit Hwang 
% last update: 06/20/2016

% Usage
%   [G,adjacency_matrix] = find_graph(edge)
%   plot(G) % you can plot a constructed graph
%   
% input: vertex - 
%        pseudoknots - 0 for no pseudoknots 
%                      1 for pseudoknots
% output: edge - 
%


% Constructs a graph with edges specified by the node pairs s and T
	x = cell2mat(edge);
	s = x(:,1);
	t = x(:,2);
	G = graph(s,t);
	A = adjacency(G);
	adjacency_matrix = full(A);

