#include <algorithm>
#include <iostream>
#include <map>
#include <stack>
#include <sstream>
#include <climits>
#include <vector>
#include <set>
#include <deque>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include <signal.h>
#include <unistd.h>
#include <cstring>

#define LEVEL 0 // 1 - print all 0 -- submission level
#define DEBUG if(LEVEL) 
volatile sig_atomic_t tle = 0;

using namespace std;

class Edge { 
	
public:
	int to;
	int length; 
	
	Edge(){}
	~Edge(){}
	Edge(int t, int l){
		to = t; length = l;
	}
	bool operator < (const Edge& e){
		return length < e.length;
	}
};

vector< vector<Edge> > minG;
int minCost=INT_MAX;

void printEdgeList(const vector< vector<Edge> > &graph, bool withWeight=false, bool isViz= false){
	for(int i=0, endI = graph.size(); i < endI; i++){
		for(int j=0, endJ = graph[i].size(); j < endJ; j++){
			if(i < graph[i][j].to){
				//~ cout << i << " -- "<< e.to << ": " << e.length << endl;
				cout << i << " "<< graph[i][j].to <<  endl;
				//~ if(withWeight){
					//~ cout << i << " "<< graph[i][j].to << " : " << graph[i][j].length <<  endl;
				//~ }else if(isViz){
					//~ cout << i << " -- "<< graph[i][j].to << "[label=" << graph[i][j].length << ",weight="<<  graph[i][j].length << ",color=red, penwidth=2]" <<  endl;
				//~ }else {
					//~ cout << i << " "<< graph[i][j].to <<  endl;
				//~ }
			}
		}
	}
	
}


int getGraphWeight(const vector< vector<Edge> > &graph){
	int mstVal =0;
	
	for(int i=0, endI = graph.size(); i < endI; i++){
		for(int j=0, endJ = graph[i].size(); j < endJ; j++){
			if(i < graph[i][j].to)
				mstVal += graph[i][j].length;
		}
	}
	
	return mstVal;
}

void term(int signum)
{
	cout << "VALUE " << minCost << "\n";
	printEdgeList(minG, false,false);
	exit(0);
} 

void addTreeEdges(set<pair<int,int>> &allEdges, const set<pair<int,int>> &edgeSet, set<int>& stopSet){
	for(auto edge : edgeSet){
		int u  = edge.first;
		int v  = edge.second;
		allEdges.insert(edge);
	
		stopSet.insert(u); 
		stopSet.insert(v); 

	}
}
void addTreeEdgesDeleteDq(set<pair<int,int>> &allEdges, 
const set<pair<int,int>> &edgeSet, 
set<int>& stopSet,
deque<int> &activeDq
){
	for(auto edge : edgeSet){
		int u  = edge.first;
		int v  = edge.second;
		allEdges.insert(edge);
	
		stopSet.insert(u); 
		stopSet.insert(v); 
		
		activeDq.erase(remove(activeDq.begin(), activeDq.end(), u), activeDq.end());
		activeDq.erase(remove(activeDq.begin(), activeDq.end(), v), activeDq.end());

	}
}
void populateEdgesOf(int u, const vector<int> & parent,  set<pair<int,int>> &edgeSet){
	int v = parent[u];			
	//~ DEBUG cout << "FOUND "<< u<<endl;
	while( v != -1){
		//~ DEBUG cout << "\t"<<  v << endl;
		if(u < v) // for speedup!
			edgeSet.insert({u,v});
		else
			edgeSet.insert({v,u});
		
		u = v;
		v =  parent[u];
	}
}
set<pair<int,int>>  getEdgesOf(int u, const vector<int> & parent){
	int v = parent[u];			
	set<pair<int,int>> edgeSet;
	//~ DEBUG cout << "FOUND "<< u<<endl;
	while( v != -1){
		//~ DEBUG cout << "\t"<<  v << endl;
		if(u < v) // for speedup!
			edgeSet.insert({u,v});
		else
			edgeSet.insert({v,u});
		
		u = v;
		v =  parent[u];
	}
	return edgeSet;
	 
}
pair<int , set<pair<int,int>> > dijkstra(
	const vector< vector<Edge> > &graph, 
	int source, 
	vector<int>& min_distance, 
	set <int> & stopSet,
	bool FirstTime = false
	) {
	
	//~ DEBUG cout << "in DJ : " << source << endl;
	set<pair<int,int>> edgeSet;
	
	vector<int> parent(graph.size() , -1);
	
	min_distance[ source ] = 0;
	set< pair<int,int> > active_vertices;
	active_vertices.insert( {0,source} );
	int where;
	set <int> nStopSet;
	
	while (!active_vertices.empty() ) {
		where = active_vertices.begin()->second;
		//~ DEBUG cout << where<<endl;
		if(stopSet.find(where) != stopSet.end()){
			
			int u = where; 
			int v = parent[u];
			nStopSet.insert(u);
			
			//~ DEBUG cout << "FOUND "<< where<<endl;
			while( v != -1){
				//~ DEBUG cout << "\t"<<  v << endl;
				nStopSet.insert(v);
				if(u < v) // for speedup!
					edgeSet.insert({u,v});
				else
					edgeSet.insert({v,u});
				
				u = v;
				v =  parent[u];
			}
			break;
		}
		
		active_vertices.erase( active_vertices.begin() );
		//cout << where<< endl ;
		for (auto ed : graph[where]) {
			auto newdist = min_distance[where] + ed.length;
			if (newdist < min_distance[ed.to]) {
				active_vertices.erase( { min_distance[ed.to], ed.to } );
				min_distance[ed.to] = newdist;
				parent[ed.to] = where;
				active_vertices.insert( { newdist, ed.to } );
			}
		}
	}
	
	//~ DEBUG cout << "end DJ : " << where << endl;
	if(FirstTime){
		stopSet.clear(); 
		stopSet.insert(nStopSet.begin(), nStopSet.end());
		//~ for (auto w : nStopSet)
			//~ cout << "== " << w<< endl;
	}else{	
		stopSet.insert(nStopSet.begin(), nStopSet.end());
	}
	return {where,edgeSet};

}
/***
 * ARG minDis is pass by ref
 * 
 * RETURNS
 * pair 1
 * =======
 * return first vertex[stopVertex] in StopSet
 * parent array
 * pair 2
 * =======
 * nonTerminals seen till then
 * Edges set for path src to the stopVertex
 */ 
//~ pair< pair<int,vector<int>> , pair< set<int>,set<pair<int,int>>> > dijkstra(
pair< pair<int,vector<int>> , pair< set<int>,set<pair<int,int>>> > dijkstra(
	const vector< vector<Edge> > &graph, 
	const set<int>& terminalSet,
	int source, 
	vector<int>& min_distance, 
	const set <int> & stopSet	
	) {
	set<int> nTermSet;
	//~ DEBUG cout << "in DJ : " << source << endl;
	set<pair<int,int>> edgeSet;
	
	vector<int> parent(graph.size() , -1);
	
	min_distance[ source ] = 0;
	set< pair<int,int> > active_vertices;
	active_vertices.insert( {0,source} );
	int where;
	
	while (!active_vertices.empty() ) {
		where = active_vertices.begin()->second;
		//~ DEBUG cout << where<<endl;
		
		if(stopSet.find(where) != stopSet.end()){
			// for optimization
			//~ populateEdgesOf(where, parent, edgeSet);
			break;
		}
		/// we are interested in the best steiner vertex
		/// That could be a terminal also!
		//~ if(terminalSet.find(where) == terminalSet.end()){
			nTermSet.insert(where);
		//~ }
		
		active_vertices.erase( active_vertices.begin() );
		//cout << where<< endl ;
		for (auto ed : graph[where]) {
			auto newdist = min_distance[where] + ed.length;
			if (newdist < min_distance[ed.to]) {
				active_vertices.erase( { min_distance[ed.to], ed.to } );
				min_distance[ed.to] = newdist;
				parent[ed.to] = where;
				active_vertices.insert( { newdist, ed.to } );
			}
		}
	}
	
	//~ DEBUG cout << "end DJ : " << where << endl;
	return {{where,parent},{nTermSet, edgeSet}};

}

/***
 * RETURNS
 * cost of the steiner tree
 * Graph/Tree DS
 */
pair<int,vector< vector<Edge> >> heuristic_snake(
	vector< vector<Edge> > &graph, 
	map<pair<int,int> , int>& W,
	set<int>& terminalSet,
	deque<int> activeDq // this has some permutation of the Terminals!
	){

	int N = graph.size();
	vector< vector<Edge> > nG(N);
	set<pair<int,int>> allEdges;
		
	set<int> stopSet;  
	 
	
	int v = activeDq.front(); activeDq.pop_front(); // first terminal
	
	stopSet.insert(v);	
	//~ DEBUG { cout << "DJ " ; printS(stopSet); }

	while(activeDq.size() > 1 ) { // take two pair at time
		int v1 = activeDq.front(); activeDq.pop_front();
		int v2 = activeDq.front(); activeDq.pop_front();
		int aCost = 0; // aCost= minD1[v1_] + minD2[v2_] ;

		vector<int> minD1(N, INT_MAX);
		auto p1 = dijkstra(graph, terminalSet, v1, minD1, stopSet);
		set<int> nonTermSet1 = p1.second.first;
		int v1_ = p1.first.first;	// read as v1 dash or v1'
		auto E1 = getEdgesOf(v1_, p1.first.second);
		auto parent1 = 	p1.first.second;
		
		vector<int> minD2(N, INT_MAX); 
		auto p2 = dijkstra(graph, terminalSet, v2, minD2, stopSet);
		int v2_ = p2.first.first; // read as v2 dash or v2'
		auto parent2 = p2.first.second;
		auto E2 = getEdgesOf(v2_, p2.first.second);				

		aCost= minD1[v1_] + minD2[v2_];
		
		set<int> nonTermSet2 = p2.second.first;
		set<int> nonTermSet;  
		//~ FIND INTESECTION
		//~ nonTermSet = nonTermSet1 \cap nonTermSet2
		set_intersection( nonTermSet1.begin(),nonTermSet1.end(),
						  nonTermSet2.begin(),nonTermSet2.end(),
					  std::inserter(nonTermSet,nonTermSet.begin()));
					 
		//~ DEBUG printS(nonTermSet);    		
		
		int minC = INT_MAX;
		int minU = -1;
		set<pair<int,int>> minE;
		
		if(nonTermSet.empty()){
			addTreeEdgesDeleteDq(allEdges, E1, stopSet, activeDq);
			addTreeEdgesDeleteDq(allEdges, E2, stopSet, activeDq);		
			continue;
		}
			
		// this is else part. i.e INTERSECTION is not empty
		for(auto u : nonTermSet ){
			vector<int> minD(N, INT_MAX);
			auto p = dijkstra(graph, terminalSet, u, minD, stopSet);
			int u_ = p.first.first; // u' is shortest from u to T
			int cost = minD[u_] + minD1[u] + minD2[u];
			auto E = getEdgesOf(u_, p.first.second);
			if(cost < minC){
				minC = cost;
				minU = u;
				minE = E;
			}
		}
		if(aCost < minC){ // claw is better for T
			addTreeEdgesDeleteDq(allEdges, E1, stopSet, activeDq);
			addTreeEdgesDeleteDq(allEdges, E2, stopSet, activeDq);
			
		}else{ //2 K2 is better for T
			
			set<pair<int,int>>edges;
			
			populateEdgesOf(minU,parent1,edges);
			populateEdgesOf(minU,parent2,edges);
			
			addTreeEdgesDeleteDq(allEdges, edges, stopSet, activeDq);	
			addTreeEdgesDeleteDq(allEdges, minE, stopSet, activeDq);
			
		}
    }
    
    if(activeDq.size()>0){ // final terminal left!
		vector<int> minD1(N, INT_MAX);
		int v1 = activeDq.front(); activeDq.pop_front();
		
		auto p1 = dijkstra(graph, terminalSet, v1, minD1, stopSet);
		auto u = p1.first.first;

		auto E1 = getEdgesOf(u, p1.first.second);
		
		addTreeEdgesDeleteDq(allEdges, E1, stopSet, activeDq);
		
	}          
	
	int mstVal = 0;
	for(auto a : allEdges){
		int u = a.first;
		int v = a.second;
		int w = W[{u,v}];
		nG[u].push_back(Edge(v, w));
		nG[v].push_back(Edge(u, w));
		mstVal += w;
		
	}
	
	//~ DEBUG cout << "MST VAL IN FUN " << mstVal << endl;
	return {mstVal,nG};
	
}

/*
 * 	RETURNS ONLY THE V1 and V2
 *  if there is no V2 then -1 
 */
pair<int, int> dijkstraSet(
	const vector< vector<Edge> > &graph, 
	const set<int>& terminalSet,
	set<int> sourceSet,  
	set <int> & stopSet
	) {
	set<int> nTermSet;
	int N = graph.size();
	
	vector<int> min_distance(N, INT_MAX);
	vector<int> parent(N , -1);
	set< pair<int,int> > active_vertices;
	for(auto source : sourceSet) {
		min_distance[ source ] = 0;
		active_vertices.insert( {0,source} );
	}
	//~ DEBUG {
		//~ cout << "in DJ : "; printS(sourceSet) ;
	//~ }
	int where=-1, whereLast=-1;
	int countStop = 0; 	int VERTEX_LIMIT = 2;
	while (!active_vertices.empty() && !stopSet.empty() ) {
		where = active_vertices.begin()->second;
		//~ DEBUG cout << where<<endl;
		if(stopSet.find(where) != stopSet.end() ){
			// for optimization
			//~ populateEdgesOf(where, parent, edgeSet);
			countStop++;
			
			stopSet.erase(where);
			//~ cout << "Removing:" << where << endl;
			if(countStop >= VERTEX_LIMIT){
				return {whereLast,where};
			}
			else if(stopSet.empty()){
				return {whereLast,-1};
			}
			whereLast = where;
		}
		
		/// we are interested in the best steiner vertex
		/// That could be a terminal also!
		//~ if(terminalSet.find(where) == terminalSet.end()){
			//~ nTermSet.insert(where);
		//~ }
		
		active_vertices.erase( active_vertices.begin() );
		//cout << where<< endl ;
		for (auto ed : graph[where]) {
			auto newdist = min_distance[where] + ed.length;
			if (newdist < min_distance[ed.to]) {
				active_vertices.erase( { min_distance[ed.to], ed.to } );
				min_distance[ed.to] = newdist;
				parent[ed.to] = where;
				active_vertices.insert( { newdist, ed.to } );
			}
		}
	}
	
	//~ DEBUG cout << "end DJ : " << where << endl;
	
	return {-1,-1};

}
vector< vector<Edge> > heuristic(
	vector< vector<Edge> > &graph, 
	map<pair<int,int> , int>& W,
	set<int>& terminalSet,
	int v  
	){
	//~ DEBUG cout << "In H1 "<< terminalSet.size() << endl;	
	int N = graph.size();
	vector< vector<Edge> > nG(N);
	set<pair<int,int>> allEdges;
	
	set<int> stopSet(terminalSet.begin(), terminalSet.end());  
	set<int> active(terminalSet.begin(), terminalSet.end());
	
	stopSet.erase(v);
	//~ int v = *(active.begin());
	
	vector<int> min_distance(N, INT_MAX);
	//~ DEBUG cout << "DJ" << endl;
	auto paired = dijkstra(graph, v,   min_distance, stopSet, true);
	int u = paired.first;
	set<pair<int,int>>  edgeSet = paired.second;
	allEdges.insert(edgeSet.begin(), edgeSet.end());
	
	
	active.erase(v);
	active.erase(u);
	//~ DEBUG cout << "DJ - done" << endl;
	while(!active.empty()){
		int v = *(active.begin());
		active.erase(active.begin());
		vector<int> min_distance(N, INT_MAX);
		
		auto paired = dijkstra(graph, v, min_distance, stopSet);
		
		auto edgeSet = paired.second;
		allEdges.insert(edgeSet.begin(), edgeSet.end());
		
	}
	
	int mstVal = 0;
	for(auto a : allEdges){
		int u = a.first;
		int v = a.second;
		int w = W[{u,v}];
		nG[u].push_back(Edge(v, w));
		nG[v].push_back(Edge(u, w));
		mstVal += w;
	}
	//~ DEBUG cout << "MST VAL " << mstVal << endl;
	return nG;
	
}

/***
 * RETURNS
 * cost of the steiner tree
 * Graph/Tree DS
 */
pair<int,vector< vector<Edge> >> nsn_fpt(
	vector< vector<Edge> > &graph, 
	map<pair<int,int> , int>& W,
	set<int>& terminalSet,
	int v // this from some terminal
	){

	int N = graph.size();
	vector< vector<Edge> > nG(N);
	set<pair<int,int>> allEdges;
	
	set<int> remTerm(terminalSet.begin(), terminalSet.end());
	remTerm.erase(v); // removing the first v!
	
	set<int> solSet;
	set<int> stopSet;
	solSet.insert(v); // solution set or src set!
	stopSet.insert(v);
	
	while(remTerm.size() > 1){
		auto dPair = dijkstraSet(graph,terminalSet, solSet, remTerm);
		int v1 = dPair.first; 
		int v2 = dPair.second;
		//~ cout << "v1:" << v1 << " v2:" << v2 << endl;
		int aCost = 0; // aCost= minD1[v1_] + minD2[v2_] ;

		vector<int> minD1(N, INT_MAX);
		auto p1 = dijkstra(graph, terminalSet, v1, minD1, stopSet);
		set<int> nonTermSet1 = p1.second.first;
		int v1_ = p1.first.first;	// read as v1 dash or v1'
		auto E1 = getEdgesOf(v1_, p1.first.second);
		auto parent1 = 	p1.first.second;
		
		vector<int> minD2(N, INT_MAX); 
		auto p2 = dijkstra(graph, terminalSet, v2, minD2, stopSet);
		int v2_ = p2.first.first; // read as v2 dash or v2'
		auto parent2 = p2.first.second;
		auto E2 = getEdgesOf(v2_, p2.first.second);				

		aCost= minD1[v1_] + minD2[v2_];
		
		set<int> nonTermSet2 = p2.second.first;
		set<int> nonTermSet;  
		//~ FIND INTESECTION
		//~ nonTermSet = nonTermSet1 \cap nonTermSet2
		set_intersection( nonTermSet1.begin(),nonTermSet1.end(),
						  nonTermSet2.begin(),nonTermSet2.end(),
					  std::inserter(nonTermSet,nonTermSet.begin()));
					 
		//~ DEBUG printS(nonTermSet);    		
		
		int minC = INT_MAX;
		int minU = -1;
		set<pair<int,int>> minE;
		
		if(nonTermSet.empty()){
			addTreeEdges(allEdges, E1, stopSet);
			addTreeEdges(allEdges, E2, stopSet);		
			continue;
		}
			
		// this is else part. i.e INTERSECTION is not empty
		for(auto u : nonTermSet ){
			vector<int> minD(N, INT_MAX);
			auto p = dijkstra(graph, terminalSet, u, minD, stopSet);
			int u_ = p.first.first; // u' is shortest from u to T
			int cost = minD[u_] + minD1[u] + minD2[u];
			auto E = getEdgesOf(u_, p.first.second);
			if(cost < minC){
				minC = cost;
				minU = u;
				minE = E;
			}
		}
		if(aCost < minC){ // claw is better for T
			addTreeEdges(allEdges, E1, stopSet );
			addTreeEdges(allEdges, E2, stopSet );
			
		}else{ //2 K2 is better for T
			
			set<pair<int,int>>edges;
			
			populateEdgesOf(minU,parent1,edges);
			populateEdgesOf(minU,parent2,edges);
			
			addTreeEdges(allEdges, edges, stopSet );	
			addTreeEdges(allEdges, minE, stopSet );
			
		}
    }
    
    if(remTerm.size()>0){ // final terminal left!
		vector<int> minD1(N, INT_MAX);
		int v1 = *(remTerm.begin()); 
		
		remTerm.erase(remTerm.begin());
		
		auto p1 = dijkstra(graph, terminalSet, v1, minD1, stopSet);
		auto u = p1.first.first;

		auto E1 = getEdgesOf(u, p1.first.second);
		
		addTreeEdges(allEdges, E1, stopSet);
		
	}          
	
	int mstVal = 0;
	for(auto a : allEdges){
		int u = a.first;
		int v = a.second;
		int w = W[{u,v}];
		nG[u].push_back(Edge(v, w));
		nG[v].push_back(Edge(u, w));
		mstVal += w;
		
	}
	
	//~ DEBUG cout << "MST VAL IN FUN " << mstVal << endl;
	return {mstVal,nG};
	
}

int main(){
	struct sigaction action;
    memset(&action, 0, sizeof(struct sigaction));
    action.sa_handler = term;
    sigaction(SIGTERM, &action, NULL);
 
    ios_base::sync_with_stdio(false);
    
	vector< vector<Edge> > graph;
	map<pair<int,int> , int> W;
	vector <int> terminals;
	set <int> terminalSet;
	deque<int> terminalsDq;
	string code, type, dummy;
	int N;
	 
	while( cin>> code >> type ){
		//cout << code  << endl;
		
		if(code == "SECTION" && type =="Graph"){
			long m, n;
			long u, v, w;
			cin >> dummy >> n;
			cin >> dummy >> m;
			N = n + 1;
			terminals.resize(N);
			
			
			//cout <<"n="<< n <<";"<< "m="<<m << endl;
			//graph= new Graph(n,true);
			graph.resize(n+1); // coz graph has from index 0. where as challege its 1
			for(long i=0; i < m; i++){
				cin>> dummy >> u >> v >> w;
				graph[u].push_back(Edge(v,w));
				graph[v].push_back(Edge(u,w));
				W[make_pair(u,v)]=w;
				W[make_pair(v,u)]=w;
				//cout << u<< " -- "<< v << " :"<< w << endl;
			}
			cin >> dummy;
		}
		else if(code == "SECTION" && type =="Terminals"){
			long t, u;
			cin >> dummy >> t;
			for(long i=0; i < t; i++){
				cin>> dummy >> u;
				//cout << "T" << u << endl;
				terminals.push_back(u);
				terminalsDq.push_back(u);
				terminalSet.insert(terminalSet.end(), u);
			}
			cin >> dummy;
		}
		else if(code == "SECTION" && type =="Tree"){
			
			cin >> dummy >> dummy >> dummy;
			long b, val ; cin >> b; 
			
			cin >> dummy >> dummy >> ws;
			
			for(long i=0; i < b; i++){
				string line;
				getline(cin, line); stringstream sstream(line);
				if(sstream >> dummy, dummy=="b"){
					while(sstream >> val){
						//cout << val << " " ;
					}
					//cout << endl;
				}
			}
			long tu, tv;
			for(long i=0; i < b-1; i++){ // b-1 edges is Td
				cin >>  tu >> tv;
				//cout<<  tu << " " << tv << endl;
			}
			cin >> dummy; // end
		}
		else{
			//cout << "Err" << endl;
		}
	
	}
	minG.resize(graph.size());
	
	/// RUN 
	/// STICKY TONGUE From all terminals
	for(auto v : terminalSet) {
		auto nG= heuristic(graph, W,terminalSet, v);
		int cost = getGraphWeight(nG);
		DEBUG cout << "COST " << minCost<< endl;
		if(cost < minCost){
			minCost = cost;
			minG = nG;	
			
		}
		
	}
	
	
	/// RUN 
	/// closest two terminals from each terminals
	
	for(auto v : terminalSet){
		auto nGPair= nsn_fpt(graph, W,terminalSet, v); // terminalsDq[0]
		int cost = nGPair.first;
		if(cost < minCost){
			minCost=cost;
			minG = nGPair.second;
		}
	}
	
	cout << "VALUE " << minCost << "\n";
	printEdgeList(minG, false,false);
	
	return 0;
}
