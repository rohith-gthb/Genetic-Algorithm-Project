#include <iostream>
#include <fstream>
#include <string>
#include<bits/stdc++.h>
#include <cmath>
using namespace std;


double h=100,w=70,lip=20,t=1,nd=5,dt=5,L=500;
double E = 210000;
double ll=20;


double centeroid(string web_half, string flange){
    // cout<<web_half<<" "<<flange<<endl;
    double ax=0,a=0;
    // double nd=5,dt=5,h=45,w=60,t=1,ll=15;
    double x=web_half[0]=='0'? 0: nd;
    
    for(int i=0; i<web_half.size()-1; i++){
        if(web_half[i]!=web_half[i+1]) {
            // cout<<"at index "<<i<<" "<<web_half[i]<<"&"<<web_half[i+1]<<"."<<endl;
            if(web_half[i]=='0'){
                x+=nd/2;
                ax+=(sqrt(nd*nd+dt*dt)*t*x);
                a+=sqrt(nd*nd+dt*dt)*t;
                x+=nd/2;
            }
            else {
                x-=nd/2;
                ax+=(sqrt(nd*nd+dt*dt)*t*x);
                a+=sqrt(nd*nd+dt*dt)*t;
                x-=nd/2;
            }
        }
        else {
            if(x!=0) ax+=(dt*t*x);
            a+=(dt*t);
        }
    }
    double wa = a;
    // cout<<"x at end of web: "<<x<<endl;
    // cout<<"web area: "<<wa<<" ";
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
    // cout<<"x at end of flange: "<<x<<endl;
    // cout<<"flange a: "<<a-wa<<" ";
    ax=ax+ll*t*w;
    a+=ll*t;
    // cout<<"tot area: "<<a<<endl;
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
    double num = M_PI*M_PI*E*I*1e-6;
    double den = L*L;
    return (num/den);
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

pair<string,string> crossover(string s, string ss){
    int n=s.size(), cr=rand()%n;
    string pre1="", pre2="", suf1="", suf2="";
    pre1=s.substr(0,cr+1); pre2=ss.substr(0,cr+1);
    suf1=s.substr(cr+1,n); suf2=ss.substr(cr+1,n);
    return {pre1+suf2, pre2+suf1}; 
}

unordered_map<string,int> dp;
int active_nodes(string s) {
    if(dp[s]!=0) return dp[s];
    int c=0;
    for(char t: s) if(t=='1') c++;
    return dp[s]=c;
}

void write_csv(string filename, vector<pair<double,pair<string,string>>> dataset){
    // Create an output filestream object
    ofstream myFile(filename);
    
    // Send column names to the stream
    myFile << "Active nodes"<<","<<"Buckling load MN";
    myFile << "\n";
    
    // Send data to the stream
    for(int i=0; i<dataset.size(); i++) {
        myFile<<active_nodes(dataset[i].second.first)+active_nodes(dataset[i].second.second)<<",";
        myFile<<dataset[i].first<<"\n";
    }
    
    myFile.close();
}

void write_best_to_text(string filename, unordered_map<int,pair<double,pair<string,string>>> map){
    // Create an output filestream object
    ofstream myFile(filename+".csv");
    myFile << "Active nodes"<<","<<"Buckling load MN";
    myFile << "\n";
    // Send data to the stream
    for(auto it: map) {
        int active_nodes=it.first;
        double buck_load=it.second.first;
        myFile<<active_nodes<<", ";
        myFile<<buck_load<<"\n";
    }
    myFile.close();
    ofstream File(filename+".txt");
    for(auto it: map) {
        int active_nodes=it.first;
        double buck_load=it.second.first;
        string web_half = it.second.second.first;
        string flange = it.second.second.second;
        File<<active_nodes<<" ";
        File<<buck_load<<" ",
        File<<web_half<<" ",
        File<<flange<<"\n";
    }

    File.close();
}


double get_buckling_load(string temp_web_half, string temp_flange) {
    double Ixx = calc_Ixx(temp_web_half, temp_flange);
    double Iyy = calc_Iyy(temp_web_half, temp_flange);
    double I = min(Iyy, Ixx);
    double buck_load=bucklingLoad(E,I,L);
    // cout<<buck_load<<endl;
    return buck_load;
}

void GA(vector<string>&web_half_vect, vector<string>&flange_vect, int iterations, double h, double w, double lip, double t, double L, double E){

    int web_size = web_half_vect[0].size(), flange_size = flange_vect[0].size();

    int popn = web_half_vect.size();
    double global_best=0; string best_web_half, best_flange;
    
    vector<pair<double,pair<string,string>>> all_points,best_in_best;

    unordered_map<int,pair<double, pair<string,string>>> map;

    int unit = iterations/100;
    for(int iter=1; iter<=iterations; iter++){
        // if(iter%unit==0) cout<<iter/unit<<"/"<<100<<endl;
        cout<<iter<<endl;
        vector<pair<double,pair<string,string>>> personal_best;
        vector<string> new_web_half_vect; 
        vector<string> new_flange_vect;
        
        for(int i=0; i<popn; i++){
            string temp_web_half = web_half_vect[i];
            string temp_flange = flange_vect[i];
            
            double buck_load = get_buckling_load(temp_web_half, temp_flange);
            
            personal_best.push_back({buck_load, {temp_web_half, temp_flange}});

            int activenodes = active_nodes(temp_web_half)+active_nodes(temp_flange);
            if(map.find(activenodes)!=map.end()){
                if(map[activenodes].first<buck_load) {
                    map[activenodes]={buck_load, {temp_web_half, temp_flange}};
                }
            }
            else {
                map[activenodes]={buck_load, {temp_web_half, temp_flange}};
            }

        	if(global_best<buck_load) {
        	    global_best=buck_load;
        	    best_web_half=temp_web_half;
        	    best_flange = temp_flange;
                best_in_best.push_back({global_best, {best_web_half ,best_flange}});
        	}
        }
        if(iter==1) {
            write_csv("initial.csv", personal_best);
        }
        // sort(personal_best.begin(), personal_best.end());

        pair<string,string> temp_best1 = personal_best[popn-1].second;
        pair<string,string> temp_best2 = personal_best[popn-2].second;
        
        vector<pair<double, pair<string,string>>> new_pers_best;
        
        for(int i=0; i<popn/2; i++){
            //take 4 parents at random
            int p1 = rand()%popn, p2=rand()%popn;
            int p3 = rand()%popn, p4=rand()%popn;
            //select two among them after competing
            int sp1, sp2;
            if(personal_best[p1].first>personal_best[p2].first) sp1=p1;
            else sp1=p2;
            if(personal_best[p3].first>personal_best[p4].first) sp2=p3;
            else sp2=p4;
            //crossover for the selected parents. 
            pair<string, string> p  = crossover( personal_best[i].second.first,  personal_best[i+popn/2].second.first);
            pair<string, string> pp = crossover(personal_best[i].second.second, personal_best[i+popn/2].second.second);
            pair<string,string> child1 = {p.first, pp.first};
            pair<string, string> child2 = {p.second, pp.second};
            double buck_load_1 = get_buckling_load(child1.first, child1.second), buck_load_2=get_buckling_load(child2.first, child2.second);
            if(buck_load_1<min(personal_best[sp1].first, personal_best[sp2].first)) {i--; continue;}
            else if(buck_load_2<min(personal_best[sp1].first, personal_best[sp2].first)) {i--; continue;}
            new_pers_best.push_back({buck_load_1,child1});
            new_pers_best.push_back({buck_load_2,child2});
        }
        priority_queue<pair<double,pair<string,string>>, 
        vector<pair<double,pair<string,string>>>, 
        greater<pair<double,pair<string,string>>>> pq;
        for(int i=0; i<popn; i++) pq.push(personal_best[i]);
        for(int i=0; i<popn; i++) {
            pq.push(new_pers_best[i]);
            pq.pop();
        }
        
        for(int i=0; i<popn; i++) {
            pair<double,pair<string,string>> p =pq.top();
            web_half_vect[i]=p.second.first;
            flange_vect[i]=p.second.second;
            pq.pop();
        }

        //Mutation 
        if(rand50()){
            int rv = rand()%popn;
            web_half_vect[rv] = mutation(web_half_vect[rv]);
            rv=rand()%popn;
            flange_vect[rv] = mutation(flange_vect[rv]);
        }

        all_points.insert(all_points.end(), personal_best.begin(), personal_best.end());
        if(iter==iterations) {
            write_csv("Final.csv", personal_best);
        }
    }
    write_csv("intermediate.csv", all_points);
    cout<<"best web half: "<<best_web_half<<" best flange: "<<best_flange<<endl;
    cout<<"global best: "<<global_best<<"MN"<<endl;
    write_csv("Best.csv", best_in_best);
    write_best_to_text("best6", map);
}



int main() {
	int popn_size=100;
	int web_half_size=(h/2)/dt+1, flange_size=w/dt+1;
	int iterations=200;
	cout<<web_half_size<<" "<<flange_size<<endl;
	vector<string> web_half_vect, flange_vect;
	srand(time(NULL));

    // Creating population 
	for(int i=1; i<=popn_size; i++) {
	    string s = "";
	    for(int j=1; j<web_half_size-1; j++) {
	        s=s+to_string(rand75_0());
	    }
	    s=s+"00";
	    web_half_vect.push_back(s);
	    s="00";
	    for(int j=3; j<flange_size-1; j++) {
	        s=s+to_string(rand75_0());
	    }
	    s=s+"00";
	    flange_vect.push_back(s);
	}

    for(int i=1; i<=popn_size; i++) {
	    string s = "";
	    for(int j=1; j<web_half_size-1; j++) {
	        s=s+to_string(rand75_1());
	    }
	    s=s+"00";
	    web_half_vect.push_back(s);
	    s="00";
	    for(int j=3; j<flange_size-1; j++) {
	        s=s+to_string(rand75_1());
	    }
	    s=s+"00";
	    flange_vect.push_back(s);
	}

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
	
	string web_half = "10000000000"; 
	string flange ="000000000000100"; 
	double Ixx = calc_Ixx(web_half, flange);
	double Iyy = calc_Iyy(web_half, flange);
	cout<<"Ixx "<<Ixx<<" Iyy: "<<Iyy<<endl;
	double I = min(Iyy, Ixx);
	cout<<"E: "<<E<<" I: "<<I<<"\n"<<"h: "<<h<<" w:"<<w<<" ll:"<<ll<<endl;
	cout<<"Original column buck load "<<bucklingLoad(E,I,L)<<"MN"<<endl;

    web_half = "10000000000"; 
	flange ="001000000000000"; 
	 Ixx = calc_Ixx(web_half, flange);
	 Iyy = calc_Iyy(web_half, flange);
	cout<<"Ixx "<<Ixx<<" Iyy: "<<Iyy<<endl;
	 I = min(Iyy, Ixx);
	cout<<"E: "<<E<<" I: "<<I<<"\n"<<"h: "<<h<<" w:"<<w<<" ll:"<<ll<<endl;
	cout<<"Original column buck load "<<bucklingLoad(E,I,L)<<"MN"<<endl;
    
    // GA(web_half_vect, flange_vect, iterations, h,w,lip,t,L,E);

	return 0;
}