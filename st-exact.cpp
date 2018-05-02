#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <algorithm>
#include <climits>
#include <signal.h>
#include <unistd.h>
#include <cstring>
#include <cmath>

#define LEVEL 0 // 1 - print all 0 -- submission level
#define DEEPDEBUG if(0) // 1 - print all 0 -- submission level
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

void printEdgeList(const vector< vector<Edge> > &graph, bool withWeight=false, bool isViz= false){
	for(int i=0, endI = graph.size(); i < endI; i++){
		for(int j=0, endJ = graph[i].size(); j < endJ; j++){
			if(i < graph[i][j].to){
			//~ cout << i << " -- "<< e.to << ": " << e.length << endl;
				if(withWeight){
					cout << i << " "<< graph[i][j].to << " : " << graph[i][j].length <<  endl;
				}else if(isViz){
					cout << i << " -- "<< graph[i][j].to << "[label=" << graph[i][j].length << ",weight="<<  graph[i][j].length  << ",weight="<<  endl;
				}else {
					cout << i << " "<< graph[i][j].to <<  endl;
				}
			}
		}
	}
	
}
template <typename T>
void printPairList(const set<pair<T,T>> & in, map<pair<int,int>,int> W, bool withWeight=false, bool isViz= false){
	for(auto p : in){ /// p - pair
		if(isViz){
			if(W.find({p.first, p.second}) != W.end()){
				int wt = W[{p.first, p.second}];
				cout << p.first << " -- " << p.second << "[label=" <<  wt  << ",color=red, penwidth=2]"<< endl;
			}
		}
		else{
			cout << p.first << " " << p.second << "\n";
		}
	}
}

// recordes the edges on path(u,v)
map<pair<int,int> , set<pair<int,int>> > dSet; 


vector <vector<int>> floydWarshall(
		vector< vector<Edge> >& graph,
		map<pair<int,int> , int> W
		){
	int n = graph.size(); // it is one plus the requried size
	vector<int> empty(n, INT_MAX/2-3); // this is ensure sum of two int-max is with in limit
	vector <vector<int>> dist(n,empty);
	 
	for(auto pairW : W){
		int u = pairW.first.first;
		int v = pairW.first.second;
		int w = pairW.second;
		//cout << u<< " " << v<< " "<< w << endl;
		dist[u][v]=w;
		dist[v][u]=w;
		pair<int,int> tmpPair; 
		if(u < v) // this is essential to ensure duplicate edge in pairs
			tmpPair = {u,v}; 
		else
			tmpPair = {v,u};
			 
		dSet[{u,v}].insert(tmpPair); // uv edge
		dSet[{v,u}].insert(tmpPair); // why?		
	}
 
	for(int k=1; k < n; k++ ){
		dist[k][k]=0;
		for(int i=1; i < n; i++ ){
			for(int j=1; j < n; j++ ){
				if(i <j && dist[i][j] > dist[i][k]+dist[k][j] ){
					dist[i][j]=dist[i][k]+dist[k][j];
					dist[j][i]=dist[i][k]+dist[k][j];
					
					dSet[{i,j}].clear(); 
					dSet[{i,j}].insert( dSet[{i,k}].begin(),  dSet[{i,k}].end() );
					dSet[{i,j}].insert( dSet[{k,j}].begin(),  dSet[{k,j}].end() );
					
					dSet[{j,i}].clear(); 
					dSet[{j,i}].insert( dSet[{i,k}].begin(),  dSet[{i,k}].end() );
					dSet[{j,i}].insert( dSet[{k,j}].begin(),  dSet[{k,j}].end() );
				}
			}
		}
	}
	
	
	return dist;
}


template <typename T>
void printV(vector<T> & in){
	cout << "[ ";
	for(auto v : in)
		cout << v << " ";
	cout << "]" << endl ;
}
template <typename T>
void printS(set<T> & in){
	cout << "( ";
	for(auto v : in)
		cout << v << " ";
	cout << ")";// << endl ;
}

template <typename T>
void printSetOfPair(set<pair<T,T>> & in){
	cout << in.size() << "=";
	for(auto pair : in)
		cout << pair.first << " " << pair.second << endl;
	
}
template <typename T>
int getWeighOfPSet(set<pair<T,T>> & in, map<pair<int,int> , int> W){
	int wt = 0;
	for(auto p : in){ /// p - pair
		wt+=W[{p.first, p.second}];
	}
	return wt;
}
map<pair<int, set<int>>,int> C; // V x 2^L -> \bb{N}
map<pair<int, set<int>>,set<pair<int,int>> > Cset; 

vector <vector<int>> d;

bool incVec(int size, vector<int>& curV){
	int cSize = curV.size(), j, max;
	for(j = cSize-1, max = size-1; j >= 0; j--) {
		if(curV[j] == max) {
			max--;
		}
		else {
			break;
		}
	}
	if(j == cSize-1) {
		curV[j]++;
	}
	else {
		while (j > 0 && curV[j] >= max) {
			j--;
			max--;
		}
		max = curV[j]+1;
		while(j < cSize) {
			curV[j++] = max++;
		}
	}
	return true;
}

int computeTab(int v, set<int> & X, vector <int> & V){
	
	int minVal =INT_MAX;
	if(X.size() ==1){
		auto x0 = *(X.begin());
		set<pair<int,int>> tmpPSet; 
		if(x0 == v){
			C.insert({{v,X}, 0});
			Cset.insert({{v,X}, tmpPSet});
			return 0;
		}	
		else{ // x0 is some u
			C.insert({{v,X}, d[v][x0]});
			tmpPSet.insert(dSet[{v,x0}].begin(), dSet[{v,x0}].end());
			
			Cset.insert({{v,X}, tmpPSet});
		}
		DEBUG {
			cout<<"T["<<v<<","; printS(X);
			cout<<"] = " <<C[{v,X}]<<"\n";
		}	
		return d[v][x0];
	}
	DEBUG {
		cout<<"computeTab("<<v<<",[";
		printS(X);
		cout<<"])\n";
	}
	vector<int> terminals(X.begin(), X.end()); // making set x as vector!
	
	set<int> X1, X2;
	X1.insert(X.begin(), X.end());
	int v1 = -1;
	int minFirst = INT_MAX;
	for(auto u: X){
		X1.erase(u);
		int cost = -1;
		if(C.find({u,X1}) == C.end()) {
			cost = computeTab(u, X1, V);
		}
		else {
			cost = C[{u,X1}];
		}
		if(minFirst > cost+d[u][v]) {
			v1 = u;
			minFirst = cost + d[u][v];
		}
		X1.insert(u);
	}
	
	int v2 = -1;
	int minSecond = INT_MAX;
	set<int> Xprime;
	int it = -1;
	for(auto u: V) {
		DEBUG{ it++; }
		if(X.find(u) != X.end()) continue;
		vector<int> curV(terminals.size(),-1); 
		int size = (1 << (terminals.size() - 1)) -1 ; // 2^(n-1) -1
		
		set<int> Xlocal;
		int minLocal = INT_MAX;
		for(int i= 0; i < size; i++) {
			incVec(terminals.size(), curV);
			X1.clear();
			for(int j = curV.size()-1; j >= 0 && curV[j] > -1; j--) {
				X1.insert(terminals[curV[j]]);
			}
			X2.clear();
			set_difference(X.begin(), X.end(), X1.begin(), X1.end(), inserter(X2, X2.begin()));
			int cost = 0;
			if(C.find({u,X1}) == C.end()) {
				cost = computeTab(u, X1, V);
			}
			else {
				cost = C[{u,X1}];
			}
			if(C.find({u,X2}) == C.end()) {
				cost += computeTab(u, X2, V);
			}
			else {
				cost += C[{u,X2}];
			}
			if(cost < minLocal) {
				minLocal = cost;
				DEBUG {
					cout<<u<<",";
					printS(X1);
					cout<<" , ";
					printS(X2);
					cout<<" , "<<cost<<endl;
				}
				Xlocal.clear();
				Xlocal.insert(X1.begin(), X1.end());
			}
 		}
 		if(minSecond > minLocal + d[u][v]) {
			minSecond = minLocal + d[u][v];
			v2 = u;
			Xprime.clear();
			Xprime.insert(Xlocal.begin(), Xlocal.end());
		}
	}
	set<pair<int, int>> P;
	set<int> Xlocal(X.begin(), X.end());
	DEBUG cout<<minFirst<<","<<minSecond<<endl;
	if(minFirst < minSecond) {
		minVal = minFirst;
		Xlocal.erase(v1);
		auto Pset = Cset[{v1,Xlocal}];
		for(auto e: Pset) {
			P.insert(e);
		}
		for(auto e: dSet[{v1,v}]) {
			P.insert(e);
		}
		DEBUG {
			cout<<"\t\tT["<<v1<<",";
			printS(Xlocal);
			cout<<"] = " << C[{v2, Xlocal}]<< " + d["<<v<<"]["<<v1<<"] = "<< d[v][v1]<<endl;
		}
	}
	else {
		minVal = minSecond;
		Xlocal.clear();
		set_difference(X.begin(), X.end(), Xprime.begin(), Xprime.end(), inserter(Xlocal, Xlocal.begin()));		
		auto Pset = Cset[{v2,Xprime}];
		for(auto e: Pset) {
			P.insert(e);
		}
		Pset = Cset[{v2,Xlocal}];
		for(auto e: Pset) {
			P.insert(e);
		}
		for(auto e: dSet[{v2,v}]) {
			P.insert(e);
		}
		DEBUG {
			cout<<"\t\tT["<<v2<<",";
			printS(Xprime);
			cout<<"] = " << C[{v2, Xprime}]<< " + T["<<v2<<",";
			printS(Xlocal);
			cout<<"] = " << C[{v2, Xlocal}] << " + d["<<v<<"]["<<v2<<"] = "<< d[v][v2]<<endl;
		}
	}
	C[{v,X}] = minVal;
	Cset.insert({{v,X}, P});
	DEBUG {
		cout<<"T["<<v<<","; printS(X);cout<<"] = " <<minVal<<"\n";
	}	
	return minVal;
}

int minVal = INT_MAX;
set<pair<int,int>> solPSet;
void term(int signum)
{
	cout << "VALUE " << minVal << endl;
	for(auto u: solPSet) {
		cout<<u.first<<" "<<u.second<< endl;
	}
	exit(0);
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
	string code, type, dummy;
	

	set<long> weightSet;
	vector<long> incident;
	 
	while( cin>> code >> type ){
		//cout << code  << endl;
		
		if(code == "SECTION" && type =="Graph"){
			long m, n;
			long u, v, w;
			cin >> dummy >> n;
			cin >> dummy >> m;
			
			incident.resize(n+1, 0) ;
			//cout <<"n="<< n <<";"<< "m="<<m << endl;
			//graph= new Graph(n,true);
			graph.resize(n+1); // coz graph has from index 0. where as challege its 1
			for(long i=0; i < m; i++){
				cin>> dummy >> u >> v >> w;
				
				weightSet.insert(w);
				incident[u]+=1;
				incident[v]+=1;
				
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
			}
			cin >> dummy; // end
		}
	}
	
	DEBUG { cout<<"graph Read\n"; }
	d = floydWarshall(graph, W); // d is global
	DEBUG { cout<<"FloydWarshall Done\n"; }
	int n = graph.size();

	set <int> terminalSet(terminals.begin(), terminals.end());
	vector <int> V;
	for(int i=1; i < n; i++)
		V.push_back(i);

	for(auto v : terminals){
		terminalSet.erase(v);
		int val = computeTab(v, terminalSet, V);
		//~ cout<<"T["<<v<<",";
		//~ printS(terminalSet);
		//~ cout<<"] = "<<val<<endl;;
		if(val < minVal){
			minVal = val;
			solPSet.clear();
			solPSet.insert(Cset[{v,terminalSet}].begin(), Cset[{v,terminalSet}].end());
		}
		terminalSet.insert(v);
	}
	cout << "VALUE " << minVal << endl;
	int y = 0;
	for(auto u: solPSet) {
		cout<<u.first<<" "<<u.second<< endl;
		y += W[{u.first,u.second}];
	}
	//~ cout << getWeighOfPSet(solPSet,W) << endl;
	//~ cout<<y<<endl;
//	printPairList(solPSet, W, false,true);
	
	return 0;
}
	
