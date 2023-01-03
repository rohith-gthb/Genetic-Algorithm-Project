#include <iostream>
#include<bits/stdc++.h>
using namespace std;

double nd=5,dt=5,h=45,w=60,t=1,ll=15;
double centeroid(string web_half, string flange){
    double ax=0,a=0;
    double nd=5,dt=5,h=45,w=60,t=1,ll=15;
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
    double dt=5,nd=5,ll=15,t=1;
    double y=web_half[0]=='0'?0:nd;
    double Ixx=0;
    double l = pow(nd*nd+dt*dt,0.5);
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
    // cout<<"y before "<<y<<endl;
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
    Ixx+=(t*pow(ll,3)/12+t*ll*pow(y-ll/2,2));
    return 2*Ixx;
}

double calc_Iyy(string web_half, string flange){
    double Xc=centeroid(web_half, flange);
    // cout<<"Xc: "<<Xc<<endl;
    double Iyy=0;
    double x=0;
    double l=pow(dt*dt+nd*nd,0.5);
    for(int i=0; i<web_half.size()-1;i++){
        if(web_half[i]!=web_half[i+1]){
            if(web_half[i]=='0'){
                x+=nd/2;
                Iyy+=(t*pow(l,3)/24+l*t*pow(Xc-x,2));
                x+=nd/2;
            }
            else{
                x-=nd/2;
                Iyy+=(t*pow(l,3)/24+l*t*pow(Xc-x,2));
                x-=nd/2;
            }
        }
        else {
            Iyy+=(pow(t,3)*dt/12+t*dt*pow(Xc-x,2));
        }
    }
    double IyyWeb=Iyy;
    for(int i=0; i<flange.size()-1;i++){
        x-=dt/2;
        if(flange[i]!=flange[i+1]){
            Iyy+=(t*pow(l,3)/24+l*t*pow(Xc-x,2));
        }
        else {
            Iyy+=(pow(t,3)*dt/12+t*dt*pow(Xc-x,2));
        }
        x-=dt/2;
    }
    double IyyFla = Iyy-IyyWeb;
    // cout<<"Iyy_web: "<<IyyWeb<<" Iyy_Flange "<<IyyFla<<" Iyy_Lip "<<(pow(t,3)*ll/12+ll*t*pow(Xc-x,2))<<endl;
    Iyy+=(pow(t,3)*ll/12+ll*t*pow(Xc-x,2));
    
    return 2*Iyy;
}

double bucklingLoad(double E, double I, double L){
    return (pow(22/7,2)*E*I/pow(L,2))/1e8;
}

int rand50()
{
	return rand() & 1;
}

string mutation(string s) {
    // cout<<"m: "<<s<<endl;
    int n = s.size();
    int i = rand()%n;
    if(s[i]=='1') s[i]='0';
    else s[i]='1';
    // cout<<"m: "<<s<<endl;
    return s;
}

pair<string,string> crossover(string s, string ss){
    // This is a single-point crossover
    // cout<<s<<ends<<ss<<endl;
    int n = s.size()/2;
    int err = s.size()%2;
    string s_new_1 = s.substr(0,n)+ss.substr(n,n+err);
    string s_new_2 = ss.substr(0,n)+s.substr(n,n+err);
    // cout<<s_new_1<<ends<<s_new_2<<endl;
    return {s_new_1, s_new_2};
}

void GA(vector<string> web_half_vect, vector<string> flange_vect, int iterations, double h, double w, double lip, double t, double L, double E){
    int popn = web_half_vect.size();
    double global_best=0; string best_web_half, best_flange;
    
    // for(string s: web_half_vect) cout<<s<<" "; cout<<endl;
    // for(string s: flange_vect) cout<<s<<" "; cout<<endl;

    for(int iter=1; iter<=iterations; iter++){

        vector<pair<double,pair<string,string>>> personal_best;
        vector<string> new_web_half_vect; 
        vector<string> new_flange_vect;
        
        for(int i=0; i<popn; i++){
            string temp_web_half = web_half_vect[i];
            string temp_flange = flange_vect[i];
            
            double Ixx = calc_Ixx(temp_web_half, temp_flange);
        	double Iyy = calc_Iyy(temp_web_half, temp_flange);
        	double I = min(Iyy, Ixx);
        	double buck_load=bucklingLoad(E,I,L);
            personal_best.push_back({buck_load, {temp_web_half, temp_flange}});
        	if(global_best<buck_load) {
        	    global_best=buck_load;
        	    best_web_half=temp_web_half;
        	    best_flange = temp_flange;
        	}
        }
        sort(personal_best.begin(), personal_best.end());
        //Crossover
        for(int i=0; i<popn/2; i++){
            pair<string, string> p = crossover(personal_best[i].second.first, personal_best[i+popn/2].second.first);
            new_web_half_vect.push_back(p.first);
            new_web_half_vect.push_back(p.second);
            pair<string, string> pp = crossover(personal_best[i].second.second, personal_best[i+popn/2].second.second);
            new_flange_vect.push_back(pp.first);
            new_flange_vect.push_back(pp.second);
        }
        //Mutation
        if(rand50()){
            mutation(new_web_half_vect[rand()%popn]);
            mutation(new_flange_vect[rand()%popn]);
        }
        
        web_half_vect=new_web_half_vect;
        flange_vect=new_flange_vect;
    }
    cout<<"best web half: "<<best_web_half<<endl;
    cout<<"best flange: "<<best_flange<<endl;
    cout<<"global best: "<<global_best<<endl;
    for(string s: web_half_vect) cout<<s<<" "<<s.size()<<endl; cout<<endl;
    for(string s: flange_vect) cout<<s<<" "<<s.size()<<endl; cout<<endl;
}

int main() {
	double h,w,lip,t,nd,dt,L;
	dt=nd=5;
	h=90; w=60; lip=15; t=1; L=75;
	double E = 210000;
	
	int popn_size=2;
	int web_half_size=10, flange_size=13;
	int iterations=10000;
	
	vector<string> web_half_vect, flange_vect;
	srand(time(NULL));
	for(int i=1; i<=popn_size; i++) {
	    string s = "0";
	    for(int j=2; j<web_half_size; j++) {
	        s=s+to_string(rand50());
	    }
	    s=s+"0";
	    web_half_vect.push_back(s);
	    s="0";
	    for(int j=2; j<flange_size; j++) {
	        s=s+to_string(rand50());
	    }
	    s=s+"0";
	    flange_vect.push_back(s);
	}
    
	// ORIGINAL COLUMN WITHOUT NOTCHES //  
	
	// string web_half = "0000000000"; 
	// string flange ="0000000000000"; 
	
	// cout<<"web_half: "<<web_half<<endl;
	// cout<<"flange: "<<flange<<endl;
	
	// double Ixx = calc_Ixx(web_half, flange);
	// double Iyy = calc_Iyy(web_half, flange);
	
	// double I = min(Iyy, Ixx);
	
	// cout<<"Original column buck load "<<bucklingLoad(E,I,L)<<"N"<<endl;
    // cout<<5*(web_half_vect[0].size()-1)<<endl;
    GA(web_half_vect, flange_vect, 1, h,w,lip,t,L,E);
    GA(web_half_vect, flange_vect, iterations, h,w,lip,t,L,E);
    
	return 0;
}