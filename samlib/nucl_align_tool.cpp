#include "nucl_align.h"
#include <string>
#include <iostream>
using namespace std;

NucleotideAlignmentTool::NucleotideAlignmentTool(){
    ;
}

NucleotideAlignmentTool::~NucleotideAlignmentTool(){
    ;
}

void NucleotideAlignmentTool::dot_matrix(NucleotideAlignment &align, vector<int> &dot, int &row, int &col){
    // restore the target sequence and query sequence by
    // remove all gaps
    string target="";
    string query="";
    int m,n;
    for (int i=0,n=0; i<(int)align.query_seq.size(); i++){
        if (align.query_seq[i]!='-'){
            query+=align.query_seq[i];
            n++;
        }
    }
    for (int i=0,m=0; i<(int)align.target_seq.size(); i++){
        if (align.target_seq[i]!='-'){
            target+=align.target_seq[i];
            m++;
        }
    }
    // resize vector
    dot.resize(m*n);
    row=m;
    col=n;
    // loop over target and query to calculate the dot matrix
    for (int i=0; i<m; i++){
        for (int j=0; j<n; j++){
            if (target[i]==query[j]){
                dot[i*n+j]=1;
            }
        }
    }
}

void NucleotideAlignmentTool::print_dot_matrix(NucleotideAlignment &align){
    vector<int> dot;
    int m,n;
    dot_matrix(align,dot,m,n);
    cout<<align.align_name<<" "<<m<<" "<<n;
    for (int i=0; i<(int)dot.size(); i++){
        cout<<" "<<dot[i];
    }
    cout<<endl;
}
