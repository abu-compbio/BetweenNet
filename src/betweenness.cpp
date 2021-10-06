// g++-5 -std=c++0x -I$LEDAROOT/incl -L$LEDAROOT betweenness.cpp -lGeoW -lD3 -lW -lP -lG -lL -lX11 -lm -O1 -no-pie -o betweenness.out

#include "common.h"
// helper method for spliting line into tokens
std::vector<std::string> splitLine(const std::string& str, const char& ch) {

	std::string next;
	std::vector<std::string> result;

	// For each character in the string
	for (std::string::const_iterator it = str.begin(); it != str.end(); it++) {
		// If we've reached the terminal character
		if (*it == ch) {
        	// If we have some characters accumulated
            	if (!next.empty()) {
                	// Add them to the result vector
                    result.push_back(next);
                    next.clear();
				}
			} else {
            	// Accumulate the next character into the sequence
                next += *it;
			}
	}
	if (!next.empty())
 		result.push_back(next);

	return result;
}

void intActDataGraph(
	string file, 									// ineraction file
	graph *intActGraph, 							// ineraction graph (to be constructed)
	std::map<int, std::string> *nodes_map_rev,		// map<index(node), gene name>
	map<int, int> *edge_weight_map,					// map<index(edge), confidence score>
	int weight_key){

	int check = 0, confidence_value;

	std::string gene, line;

	//read the intact file
	std::ifstream conv_reader(file);


	set<node> nodes_set; //set of intact nodes
	set<edge> edges_set;	//set of intact edges
	std::map<std::string, node> nodes_map;			// map<gene name, node>
	std::set<std::string> nodes_set_str; 			// all the labels
	node n, v;
	edge e;

	// getting all the labels (genes names)
	while (getline(conv_reader, line)) {
		if (check++ > 1) {
			std::vector <std::string> tokens = splitLine(line, '\t');
			nodes_set_str.insert(tokens[0]);	//add  the first interactor to the list
			nodes_set_str.insert(tokens[1]);	//add the second interactor to the list
		}
	}

	conv_reader.close();
	// creating a node for each labels
	for (auto gene : nodes_set_str)
	{
		n = (*intActGraph).new_node();
		nodes_set.insert(n);
		nodes_map[gene] = n;
		(*nodes_map_rev)[index(n)] = gene;
	}

	// creating edges
	std::ifstream conv_reader1(file);
	check = 0;
	while (getline(conv_reader1, line))
		if (check++ > 0) {
			std::vector <std::string> tokens = splitLine(line, '\t');
			n = nodes_map[tokens[0]];
			v = nodes_map[tokens[1]];

			e = (*intActGraph).new_edge(n, v);
			edges_set.insert(e);
		}
	conv_reader1.close();

	
	(*intActGraph).make_undirected();
	Make_Simple((*intActGraph));
	(*intActGraph).write_gml("..out//intact_graph.gml");
}

void construct_normal_tumor_networks(
	graph intActGraph, 								// interaction graph
	graph *normalGraph, graph *tumorGraph, 			// normal and tumor graph
	std::map<int, std::string> nodes_map_rev,		// map<index(node), node label>
	map<int, int> *normal_edge_map,					// <original edge, virtal edge>  for normal graph
	map<int, int> *tumor_edge_map,					// <original edge, virtal edge>	 for tumor graph
	std::map<int, std::string> *normal_nodes_map_rev,	// map<index(node), node label> for normal graph
	std::map<int, std::string> *tumor_nodes_map_rev,	// map<index(node), node label> for tumor graph
	string input_data,
	string filein){


	// cloning the inActGraph into *normalGraph and *tumorGraph
	edge e, f;
	node start_normal, end_normal, start_tumor, end_tumor, n, v;

	map<node, node> NormalNodeMap, TumorNodeMap;			// <virtual node, original node>
	map<node, node> NormalNodeMapRev, TumorNodeMapRev;		// <original node, virtual node>


	//create a copy of the intact ppi graph for normal network
	forall_nodes(n, intActGraph){
		v = (*normalGraph).new_node();
		NormalNodeMap[v] = n;
		NormalNodeMapRev[n] = v;
		(*normal_nodes_map_rev)[index(v)] = nodes_map_rev[index(n)];
	}

	forall_edges(e, intActGraph){
		f = (*normalGraph).new_edge(NormalNodeMapRev[intActGraph.source(e)], NormalNodeMapRev[intActGraph.target(e)]);
		(*normal_edge_map)[index(e)] = index(f);
	}

	//create a copy of the intact ppi graph for tumor network
	forall_nodes(n, intActGraph){
		v = (*tumorGraph).new_node();
		TumorNodeMap[v] = n;
		TumorNodeMapRev[n] = v;
		(*tumor_nodes_map_rev)[index(v)] = nodes_map_rev[index(n)];
	}

	forall_edges(e, intActGraph){
		f = (*tumorGraph).new_edge(TumorNodeMapRev[intActGraph.source(e)], TumorNodeMapRev[intActGraph.target(e)]);
		(*tumor_edge_map)[index(e)] = index(f);
	}


	Make_Simple((*normalGraph));
	(*normalGraph).make_undirected();
	Make_Simple((*tumorGraph));
	(*tumorGraph).make_undirected();
	// done with creating a copy of ppi netowrk
	// cloning finished


	// getting all patient normal and tumor genes
	int checkLine = 0;
	std::ifstream conv_reader(input_data + "/" + filein); //reading the file for patient p
	std::string line;
	list <std::string> gene_data_normal, gene_data_tumor; //a list to store genes in the normal graph and tumor graph
	std::vector<std::string> line_vector, words_vector;
	while (getline(conv_reader, line)) {
		if (checkLine++ > 0) {
			line_vector = splitLine(line, '\t'); //splitting the line in the file
			if (line_vector[0] != "removed") { // removed means the gene is not expressed, RSEM<5
				words_vector = splitLine(line_vector[0], '|');
				gene_data_tumor.append(words_vector[0]);
			}

			if (line_vector[2] != "removed") {
				words_vector = splitLine(line_vector[2], '|');
				gene_data_normal.append(words_vector[0]);
			}
		}
	}

	// keeping nodes belongs to normal genes list and deleting others
	cout << "Normal Graph node deletion starts\n";
	forall_nodes(n, (*normalGraph))
		if (gene_data_normal.search((*normal_nodes_map_rev)[index(n)]) == nil){ // if the node in the PPI netowrk not in normal gene set, then remove it
			(*normalGraph).del_node(n);
			(*normal_nodes_map_rev).erase(index(n));
		}
	Make_Simple((*normalGraph));
	(*normalGraph).make_undirected();
	cout << "Normal Graph # nodes " << (*normalGraph).number_of_nodes() << endl;

	// keeping nodes belongs to tumor genes list and deleting others
	cout << "Tumor Graph node deletion starts\n";
	forall_nodes(n, (*tumorGraph))
		if (gene_data_tumor.search((*tumor_nodes_map_rev)[index(n)]) == nil){
			(*tumorGraph).del_node(n);
			(*tumor_nodes_map_rev).erase(index(n));
		}
	Make_Simple((*tumorGraph));
	(*tumorGraph).make_undirected();
	cout << "Tumor Graph # nodes " << (*tumorGraph).number_of_nodes() << endl;

}

void betweenness_centrality(
	std::map<int, float> *bw,						// betweenness map<index(node), betweenness value>
	std::map<int, std::string> *nodes_map_reverse, 	// map<index(node), node label>
	map<int, int> *edge_map,						// map<original edge, virtal edge>
	graph *G 										// input graph
	) {

	// initialization
	node s, w, v, n, v_start, v_end, v_node;
	edge e, v_edge;

	std::map<int, std::string> nodes_map_rev = (*nodes_map_reverse);

	int n_nodes = (*G).number_of_nodes();

	cout << "#Nodes of (*G): " << (*G).number_of_nodes() << endl;
	cout << "#Edges of (*G): " << (*G).number_of_edges() << endl;


	cout << "Betweenness Calculation Started" << endl;
	std::map<int, float> sigma, dist, delta;
	list<int> stack;
	std::map<int, list<int>> pred;

	// initializing betweenness map to zeros for all node in the virtual graph
	forall_nodes(s, (*G))
		(*bw)[index(s)] = 0.0;

	// single shortest path calculation
	int count= 0;
	int ittt=0;
	forall_nodes(s, (*G)){
		stack.clear();
		forall_nodes(v, (*G)) {
			pred[index(v)].clear();
	 		dist[index(v)] = -1.0;		//infinity
	 		sigma[index(v)] = 0.0;
		}

		dist[index(s)] = 0.0;
 		sigma[index(s)] = 1.0;
 		list<node> queue;
		queue.append(s);
		while(!queue.empty()) {
			v = queue.pop();
			stack.append(index(v));

			float dist_v = dist[index(v)];
			forall_adj_nodes(w, v){
				if (dist[index(w)] == -1.0) {	// visited for the first time
					queue.append(w);
					dist[index(w)] = dist_v + 1.0;
				}
				if (dist[index(w)] == (dist_v + 1.0) && dist[index(w)]<=3) {
					sigma[index(w)] += sigma[index(v)];
					pred[index(w)].append(index(v));
				}
			}
		}


		// accumulation part
		int acc_w, acc_v, acc_i;
		forall(acc_i, stack){
			delta[acc_i] = 0.0;
		}

 		while(!stack.empty()){
 			acc_w = stack.pop_back();
 			float coff = (1.0 + delta[acc_w]) / sigma[acc_w];
 			forall(acc_v, pred[acc_w]){
 				delta[acc_v] += sigma[acc_v] * coff;
 			}
 			if (acc_w != index(s)){
 				(*bw)[acc_w] += delta[acc_w];
 			}
 		}
	}
	cout << "Betweenness Calculation Finished!" << endl;

}

int main(int argc, char *argv[]) {
	graph intActGraph;
	edge e;
	node n;
	string patient;


	string ppi_file = "../data/IntAct_network.txt";
	leda::string cancer_type = argv[1];
	
	string input_data = "../data/"+cancer_type+"/betweenness_input/";



	int weight_key = 1;
	int counter = 0;

	std::map<int, std::string> nodes_map_rev;
	map<int, int> edge_weight_map;

	//constructing the PPI network
	intActDataGraph(ppi_file, &intActGraph, &nodes_map_rev, &edge_weight_map, weight_key);
	cout << "\n ## 1- PPI Network Generated with " << intActGraph.number_of_nodes() << " and " << intActGraph.number_of_edges() << endl;;

	list<string> files = get_files(input_data);
	forall(patient, files) {

		graph tumorGraph, normalGraph;

		cout << "\n===============================\n";
		cout << ++counter << ".Patient " << patient << " calculations starts\n";

		map<int, int> normal_edge_map, tumor_edge_map;					// <original edge, virtal edge>

		std::map<int, std::string> normal_nodes_map_rev, tumor_nodes_map_rev;

		// constructing the normal & tumor graph
		construct_normal_tumor_networks(intActGraph,
			&normalGraph, &tumorGraph,
			nodes_map_rev,
			&normal_edge_map, &tumor_edge_map,
			&normal_nodes_map_rev, &tumor_nodes_map_rev,
			input_data,
			patient);

		cout << "Normal and Tumor Network Construction Finished\n\n";
		// calculating the betweeness

		clock_t start,end;
  	start = clock();
		std::map<int, float> bw_n, bw_t;

		//calculating betweenness: normal network
		betweenness_centrality(&bw_n, &normal_nodes_map_rev, &normal_edge_map, &normalGraph);
		cout << "Normal Betweenness Calculated\n";

		std::ofstream bw_normal_output("../out/"+cancer_type+"/Betweenness/normal_bw_" + patient);
		if (bw_normal_output.is_open())
			for (std::map<int, std::string>::iterator it = normal_nodes_map_rev.begin(); it != normal_nodes_map_rev.end(); ++it)
				if ((it -> second).find("virt") != 0)
					bw_normal_output << it -> second << "\t:\t" << bw_n[(it -> first)] << endl;

		bw_normal_output.close();

		end = clock() - start;
		cout << "Normal Graph processed in " << ((float)end)/CLOCKS_PER_SEC << " seconds.\n";


		//calculating betweenness: tumor network
  	start = clock();
		betweenness_centrality(&bw_t, &tumor_nodes_map_rev, &tumor_edge_map, &tumorGraph);
		cout << "Tumor Betweenness Calculated\n";
		std::ofstream bw_tumor_output("../out/"+cancer_type+"/Betweenness/tumor_bw_" + patient);
		if (bw_tumor_output.is_open())
			for (std::map<int, std::string>::iterator it = tumor_nodes_map_rev.begin(); it != tumor_nodes_map_rev.end(); ++it)
				if ((it -> second).find("virt") != 0)
					bw_tumor_output << it -> second << "\t:\t" << bw_t[(it -> first)] << endl;
		bw_tumor_output.close();
		end = clock() - start;
		cout << "Tumor Graph processed in " << ((float)end)/CLOCKS_PER_SEC << " seconds.\n";

	}


	return 0;
}
