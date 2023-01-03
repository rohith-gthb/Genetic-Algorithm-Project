#include <iostream>
#include <fstream>
#include <string>
#include<bits/stdc++.h>
using namespace std;

#define data pair<double, pair<string,string>>

double h=100,w=70,lip=20,t=1,nd=5,dt=5,L=500;
double E = 210000;
double ll=20;


double centeroid(string web_half, string flange){
    double ax=0,a=0;
    // double nd=5,dt=5,h=45,w=60,t=1,ll=15;
    double x=web_half[0]=='0'? 0: nd;
    
    for(int i=0; i<web_half.size()-1; i++){
        if(web_half[i]!=web_half[i+1]) {
            if(web_half[i]=='0'){
                x+=nd/2;
                ax+=(pow(nd*nd+dt*dt,1/2)*t*x);
                a+=pow(nd*nd+dt*dt,1/2)*t;
                x+=nd/2;
            }
            else {
                x-=nd/2;
                ax+=(pow(nd*nd+dt*dt,1/2)*t*x);
                a+=pow(nd*nd+dt*dt,1/2)*t;
                x-=nd/2;
            }
        }
        else {
            if(x!=0) ax+=(dt*t*x);
            a+=(dt*t);
        }
    }
    
    for(int i=0; i<flange.size()-1; i++){
        x+=dt/2;
        if(flange[i]!=flange[i+1]){
            ax+=(pow(nd*nd+dt*dt,1/2)*t*x);
            a+=(pow(nd*nd+dt*dt,1/2)*t);
        }
        else {
            ax+=(dt*t*x);
            a+=(dt*t);
        }
        x+=dt/2;
    }
    
    ax=ax+ll*t*w;
    a+=ll*t;
    
    return ax/a;
}

double calc_Ixx(string web_half, string flange){
    // double y=web_half[0]=='0'?0:nd;
    double y=0;
    double Ixx=0;
    double l = pow(nd*nd+dt*dt,0.5);
    // cout<<"exprecting web half size to be "<<(web_half.size()-1)*dt<<" with size as "<<web_half.size()<<endl;
    for(int i=0; i<web_half.size()-1; i++){
        y+=dt/2;
        if(web_half[i]!=web_half[i+1]){
            Ixx+=(t*pow(l,3)/24 + t*l*pow(y,2));
        }
        else {
            Ixx+=(t*pow(dt,3)/12+t*dt*pow(y,2));
        }
        y+=dt/2;
    }
    double IxxWeb=Ixx;
    // cout<<"y before flange: "<<y<<endl;
    for(int i=0; i<flange.size()-1; i++){
        if(flange[i]!=flange[i+1]){
            if(flange[i]=='0'){
                y-=dt/2;
                Ixx+=(t*pow(l,3)/24+t*l*pow(y,2));
                y-=dt/2;
            }
            else {
                y+=nd/2;
                Ixx+=(t*pow(l,3)/24+t*l*pow(y,2));
                y+=nd/2;
            }
        }
        else {
            Ixx+=(dt*pow(t,3)/12+t*dt*pow(y,2));
        }
    }
    double IxxFla = Ixx-IxxWeb;
    // cout<<"Ixx_web: "<<IxxWeb<<" Ixx_flange: "<<IxxFla<<" Ixx_lip "<<(t*pow(ll,3)/12+t*ll*pow(y-ll/2,2))<<endl;
    // cout<<"y at last in Ixx: "<<y<<endl;
    Ixx+=(t*pow(ll,3)/12+t*ll*pow(y-ll/2,2));
    return 2*Ixx;
}

double calc_Iyy(string web_half, string flange){
    double Xc=centeroid(web_half, flange);
    // cout<<"Xc: "<<Xc<<endl;
    double Iyy=0;
    double x=-Xc;
    double l=pow(dt*dt+nd*nd,0.5);
    // cout<<"initial x in Iyy: "<<x<<" ";
    for(int i=0; i<web_half.size()-1;i++){
        if(web_half[i]!=web_half[i+1]){
            if(web_half[i]=='0'){
                x+=nd/2;
                Iyy+=(t*pow(l,3)/24+l*t*pow(x,2));
                x+=nd/2;
            }
            else{
                x-=nd/2;
                Iyy+=(t*pow(l,3)/24+l*t*pow(x,2));
                x-=nd/2;
            }
        }
        else {
            Iyy+=(pow(t,3)*dt/12+t*dt*pow(x,2));
        }
    }
    // cout<<" at the end "<<x<<endl;
    double IyyWeb=Iyy;
    for(int i=0; i<flange.size()-1;i++){
        x+=dt/2;
        if(flange[i]!=flange[i+1]){
            Iyy+=(t*pow(l,3)/24+l*t*pow(x,2));
        }
        else {
            Iyy+=(pow(t,3)*dt/12+t*dt*pow(x,2));
        }
        x+=dt/2;
    }
    // cout<<"Before the lip x: "<<x<<endl;
    double IyyFla = Iyy-IyyWeb;
    // cout<<"Iyy_web: "<<IyyWeb<<" Iyy_Flange "<<IyyFla<<
    // " Iyy_Lip "<<(pow(t,3)*ll/12+ll*t*pow(x,2))<<endl;
    Iyy+=(pow(t,3)*ll/12+ll*t*pow(x,2));
    
    return 2*Iyy;
}

double bucklingLoad(double E, double I, double L){
    return (pow(22/7,2)*E*I/pow(L,2))/1e6;
}

int rand50()
{
	return rand() & 1;
}

int rand75_0()
{
	return rand50() & rand50();
}

int rand75_1()
{
	return 1-rand75_0();
}

string mutation(string s) {
    int n = s.size();
    int i = rand()%n;
    while(i==n-1 || i==n-2) i=rand()%n;
    if(n==15) while(i<2 || i>=n-2) i=rand()%n;
    if(s[i]=='1') s[i]='0';
    else s[i]='1';
    return s;
}

pair<string,string> crossover_single(string s, string ss){
    int n=s.size();
    int cr=rand()%n;
    string pre1="", pre2="", suf1="", suf2="";
    pre1=s.substr(0,cr+1); pre2=ss.substr(0,cr+1);
    suf1=s.substr(cr+1,n); suf2=ss.substr(cr+1,n);
    return {pre1+suf2, pre2+suf1}; 
}

void crossover(pair<string, string>&p, pair<string, string>&pp){
    pair<string,string> webs = crossover_single(p.first, pp.first);
    pair<string,string> flanges = crossover_single(p.second, pp.second);
    p.first=webs.first; pp.first=webs.second;
    p.second=flanges.first; pp.second=flanges.second;
}

// unordered_map<string,int> dp;
int active_nodes(string s) {
    // if(dp[s]!=0) return dp[s];
    int c=0;
    for(char t: s) if(t=='1') c++;
    return c;
    // return dp[s]=c;
}

void write_csv(string filename, vector<pair<int,double>> dataset){
    // Create an output filestream object
    ofstream myFile(filename+".csv");
    
    // Send column names to the stream
    myFile << "Active nodes"<<","<<"Buckling load MN";
    myFile << "\n";
    
    // Send data to the stream
    for(int i=0; i<dataset.size(); i++) {
        myFile<<dataset[i].first<<", ";
        myFile<<dataset[i].second<<"\n";
    }
    
    myFile.close();
}

void write_best_to_csv(string filename, vector<pair<double,pair<string,string>>> arr, bool only_text=false){
    
    if(!only_text) {
        ofstream myFile(filename+".csv");
        myFile << "Active nodes"<<","<<"Buckling load MN";
        myFile << "\n";
        for(auto it: arr) {
            pair<string, string> p = it.second;
            double bl = it.first;
            myFile<<active_nodes(p.first)+active_nodes(p.second)<<", "<<bl<<"\n";
        }
        myFile.close();
    }

    ofstream File(filename+".txt");
    for(auto it: arr) {
        pair<string, string> p = it.second;
        double bl = it.first;
        File<<active_nodes(p.first)+active_nodes(p.second)<<" "<<bl<<" "<<p.first<<" "<<p.second<<"\n";
    }

    File.close();

}


void sharing_function(int sigma, vector<data> &fitness) {
    vector<pair<double,double>> coords;
    for(int i=0; i<fitness.size(); i++) {
        double x = active_nodes(fitness[i].second.first)+ active_nodes(fitness[i].second.second);
        double y = fitness[i].first;
        coords.push_back({x,y});
    }
    vector<double> tot_sh(fitness.size(),0.001); // cant have this to be zero as we would be using it to divide afterwards. 
    for(int i=0; i<coords.size(); i++) {
        for(int j=i+1; j<coords.size(); j++) {
            double dij = (coords[i].first-coords[j].first)*(coords[i].first-coords[j].first)
                + (coords[i].second-coords[j].second)*(coords[i].second-coords[j].second) ;
            dij = sqrt(dij);
            if(dij>=sigma) continue;
            dij = 1-pow(dij/sigma,2);
            tot_sh[i]+=dij;
            tot_sh[j]+=dij;
        }
    }
    for(int i=0; i<fitness.size(); i++) {
        fitness[i].first/=tot_sh[i];
    }
    sort(fitness.begin(), fitness.end());
    return; 
}




void GA(vector<string>&web_half_vect, vector<string>&flange_vect, int iterations, double h, double w, double lip, double t, double L, double E){    
    double best_bl=0;
    vector<data> all;
    int denominator = iterations/100;
    // unordered_map<int,priority_queue<pair<double,int>>> map;
    for(int i=1; i<=iterations; i++) {
        if(i%denominator==0) cout<<"iter: "<< i/denominator<<"/ "<<100<<endl;
        // unordered_set<string> set;
        unordered_map<int,priority_queue<pair<double,int>>> map;
        int max_an=0, best_index=-1;
        for(int j=0; j<web_half_vect.size(); j++) {
            double I = min(calc_Ixx(web_half_vect[j], flange_vect[j]), calc_Iyy(web_half_vect[j], flange_vect[j]));
            double buck_load = bucklingLoad(E,I,L);
            int acn = active_nodes(web_half_vect[j])+active_nodes(flange_vect[j]);
            all.push_back({buck_load, {web_half_vect[j], flange_vect[j]}});
            if(max_an==0 || map[max_an].top().first<buck_load) max_an=acn;
            if(map.find(acn)!=map.end()) {
                if(map[acn].top().first==buck_load) continue;
            }
            map[acn].push({buck_load,j});
            best_bl=max(best_bl,buck_load);
            if(best_bl==buck_load) {
                best_index=j;
            }
        }
        vector<data> rank1, rank2;
        for(int j=0; j<=max_an; j++) {
            if(map.find(j)==map.end()) continue;
            int top = map[j].top().second;
            rank1.push_back({map[j].top().first,{web_half_vect[top], flange_vect[top]}});
            map[j].pop();
        }
        
        for(int j=0; j<=max_an; j++) {
            if(map.find(j)==map.end()) continue;
            if(map[j].size()==0) continue;
            int top = map[j].top().second;
            rank2.push_back({map[j].top().first,{web_half_vect[top], flange_vect[top]}});
            map[j].pop();
        }
        
        if(i==iterations) {
            // write_best_to_csv("rank1",rank1);
            // write_best_to_csv("rank2",rank2);
            // write_best_to_csv("all",all);
        }

        sharing_function(10, rank1);
        sharing_function(10, rank2);

        vector<string> new_webs, new_flanges;
        for(auto it: rank1) {
            new_webs.push_back(it.second.first);
            new_flanges.push_back(it.second.second);
        }
        for(auto it: rank2) {
            new_webs.push_back(it.second.first);
            new_flanges.push_back(it.second.second);
        }
        for(int j=0; j<100; j++) {
            int r1 = rand()%rank1.size();
            int r2 = rand()%rank1.size();
            crossover(rank1[r1].second,rank1[r2].second);
            r1 = rand()%rank2.size();
            r2 = rand()%rank2.size();
            crossover(rank2[r1].second,rank2[r2].second);
        }

        for(auto it: rank1) {
            new_webs.push_back(it.second.first);
            new_flanges.push_back(it.second.second);
        }
        for(auto it: rank2) {
            new_webs.push_back(it.second.first);
            new_flanges.push_back(it.second.second);
        }
        if(rand50()) {
            int rv = rand()%new_webs.size();
            new_webs[rv]=mutation(new_webs[rv]);
            rv = rand()%new_webs.size();
            new_flanges[rv]=mutation(new_flanges[rv]);
        }
        if(i==iterations) {
            cout<<web_half_vect[best_index]<<" "<<flange_vect[best_index]<<endl;
        }
        web_half_vect=new_webs;
        flange_vect= new_flanges;
    }
    cout<<best_bl<<"MN"<<endl;
}


void original(int whs, int fs){
    string web_half = "00000000000"; 
	string flange ="00000000"; 

    cout<<h<<" "<<w<<" "<<ll<<" "<<L<<endl;

	cout<<"web_half: "<<web_half<<endl;
	cout<<"flange: "<<flange<<endl;
	
	double Ixx = calc_Ixx(web_half, flange);
	double Iyy = calc_Iyy(web_half, flange);
	
	double I = min(Iyy, Ixx);

    // cout<<"(22/7)2 x "<<E<<"x"<<I<<"/"<<L<<"2"<<endl;
    
	cout<<"Ixx: "<<Ixx<<" Iyy: "<<Iyy<<endl;
	cout<<"Original column buck load "<<bucklingLoad(E,I,L)<<"MN"<<endl;
}


// double h=100,w=70,lip=20,t=1,nd=5,dt=5,L=500;
// double E = 210000;
// double ll=20;


int main() {
	int popn_size=50;
	int web_half_size=h/2/dt+1, flange_size=w/dt+1;
    cout<<"Given web_h_size: "<<web_half_size<<" "<<" flange_size: "<<flange_size<<endl;
	int iterations=1000;
	
	vector<string> web_half_vect, flange_vect;
	srand(time(NULL));

    // Creating population 
	for(int i=1; i<=popn_size; i++) {
	    string s = "";
	    for(int j=1; j<web_half_size-1; j++) {
	        s=s+to_string(rand50());
	    }
	    s=s+"00";
        
	    web_half_vect.push_back(s);
	    s="00";
	    for(int j=3; j<flange_size-1; j++) {
	        s=s+to_string(rand50());
	    }
	    s=s+"00";
        
	    flange_vect.push_back(s);
	}
	// ORIGINAL COLUMN WITHOUT NOTCHES //  
	original(web_half_size, flange_size);
    // Genetic algorithm
    // GA(web_half_vect, flange_vect, iterations, h,w,lip,t,L,E);
	return 0;
}