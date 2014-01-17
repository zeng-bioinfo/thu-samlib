#include "nucl_align.h"
#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "SSTree.h"
using namespace std;

NucleotideAlignmentMethod::NucleotideAlignmentMethod(){
    ;
}

NucleotideAlignmentMethod::~NucleotideAlignmentMethod(){
    ;
}

/**
 * @brief NucleotideAlignmentMethod::dot_matrix
 * @param align
 * @param dot
 * @param row
 * @param col
 */
void NucleotideAlignmentMethod::dot_matrix(NucleotideAlignment &align, vector<int> &dot, int &row, int &col){
    // restore the target sequence and query sequence by
    // remove all gaps
    string target="";
    string query="";
    int m=0,n=0;
    for (int i=0; i<(int)align.query_seq.size(); i++){
        if (align.query_seq[i]!='-'){
            query+=align.query_seq[i];
            n++;
        }
    }
    for (int i=0; i<(int)align.target_seq.size(); i++){
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

/**
 * @brief NucleotideAlignmentMethod::print_dot_matrix
 * @param align
 */
void NucleotideAlignmentMethod::print_dot_matrix(NucleotideAlignment &align){
    vector<int> dot;
    int m,n;
    dot_matrix(align,dot,m,n);

    vector<string> tokens;
    boost::split(tokens, align.align_name, boost::is_any_of(" "));

    cout<<tokens[1]<<" "<<m<<" "<<n;
    for (int i=0; i<(int)dot.size(); i++){
        cout<<" "<<dot[i];
    }
    cout<<endl;

}

/**
 * @brief NucleotideAlignmentMethod::maximal_pairs
 * @param align
 * @param p1
 * @param p2
 * @param len
 */
void NucleotideAlignmentMethod::maximal_pairs(const NucleotideAlignment &align, vector<int> &pt, vector<int> &pq, vector<int> &len){
    string raw_query=align.raw_query;
    string raw_target=align.raw_target;

    // synthesize target and query together
    char *text=new char[raw_target.size()+raw_query.size()+1+1];
    sprintf(text,"%s@%s",raw_target.c_str(),raw_query.c_str());
    ulong n=strlen(text);
    n++;    // this takes the end symbol 0u as part of the text

    // position of symbol $
    ulong splitpos=strlen(raw_target.c_str());

    // suffix tree
    SSTree *sst=new SSTree((uchar*)text,n);
    ulong lastleaf=sst->rightmost(0);


    int alpha_size=64;
    int node_size=lastleaf+1;
    // llc:left character on target
    // rlc:right character on query
    int *llc=new int[alpha_size*node_size], *rlc=new int[alpha_size*node_size];
    int *lln=new int[alpha_size*node_size], *rln=new int[alpha_size*node_size];
    int *llt=new int[alpha_size*node_size], *rlt=new int[alpha_size*node_size];

    // string_depth of nodes
    int *depth=new int[node_size];

    for (int i=lastleaf;i>0;i--) {
        depth[i]=sst->depth(i);
    }

    // finding all maximal pairs
    // bottom-up traversal
    for (int i=lastleaf; i>0; i--){
        if (sst->isleaf(i) && sst->isOpen(i)){
            // starting position of the suffix
            int ctp=sst->textpos(i);
            // left character
            int clc;
            if (ctp>0) clc=(char)text[ctp-1]-'@';
            else clc='^'-'@';
            // first string
            if (ctp<splitpos){
                llc[clc*(lastleaf+1)+i]=i;
                lln[clc*(lastleaf+1)+i]=1;
                llt[clc*(lastleaf+1)+i]=i;
            }
            // second string
            if (ctp>splitpos){
                rlc[clc*(lastleaf+1)+i]=i;
                rln[clc*(lastleaf+1)+i]=1;
                rlt[clc*(lastleaf+1)+i]=i;
            }
        }


        if (!sst->isleaf(i) && sst->isOpen(i)){
            int j=sst->firstChild(i);
            for (int x=0; x<alpha_size; x++){
                llc[x*(lastleaf+1)+i]=llc[x*(lastleaf+1)+j];
                lln[x*(lastleaf+1)+i]=lln[x*(lastleaf+1)+j];
                llt[x*(lastleaf+1)+i]=llt[x*(lastleaf+1)+j];

                rlc[x*(lastleaf+1)+i]=rlc[x*(lastleaf+1)+j];
                rln[x*(lastleaf+1)+i]=rln[x*(lastleaf+1)+j];
                llt[x*(lastleaf+1)+i]=rlt[x*(lastleaf+1)+j];
            }


            j=sst->sibling(j);
            while (j!=0){
                for (int x=0; x<alpha_size; x++){
                    for (int y=0; y<alpha_size; y++){
                        if (x==y) continue;
                        // compute the maximal pairs (left,right)
                        for (int u=0,si=llc[x*(lastleaf+1)+i]; u<lln[x*(lastleaf+1)+i]; u++,si=llc[x*(lastleaf+1)+si]){
                            int lp=sst->textpos(si);
                            for (int v=0,sj=rlc[y*(lastleaf+1)+j]; v<rln[y*(lastleaf+1)+j]; v++,sj=rlc[y*(lastleaf+1)+sj]){
                                int rp=sst->textpos(sj);
                                if (depth[i]>=5){
                                    pt.push_back(lp);
                                    pq.push_back(rp-splitpos-1);
                                    len.push_back(depth[i]);
                                }
                            }
                        }
                        // compute the maximal pairs (right,left)
                        for (int u=0,si=rlc[x*(lastleaf+1)+i]; u<rln[x*(lastleaf+1)+i]; u++,si=rlc[x*(lastleaf+1)+si]){
                            int rp=sst->textpos(si);
                            for (int v=0,sj=llc[y*(lastleaf+1)+j]; v<lln[y*(lastleaf+1)+j]; v++,sj=llc[y*(lastleaf+1)+sj]){
                                int lp=sst->textpos(sj);
                                if (depth[i]>=5){
                                    pt.push_back(lp);
                                    pq.push_back(rp-splitpos-1);
                                    len.push_back(depth[i]);
                                }
                            }
                        }
                    }
                }
                // update the list
                for (int x=0; x<alpha_size; x++){
                    if (lln[x*(lastleaf+1)+j]>0){
                        // left
                        if (lln[x*(lastleaf+1)+i]!=0){
                            llc[x*(lastleaf+1)+llt[x*(lastleaf+1)+i]]=llc[x*(lastleaf+1)+j];
                            lln[x*(lastleaf+1)+i]+=lln[x*(lastleaf+1)+j];
                            llt[x*(lastleaf+1)+i]=llt[x*(lastleaf+1)+j];
                        }else{
                            llc[x*(lastleaf+1)+i]=llc[x*(lastleaf+1)+j];
                            lln[x*(lastleaf+1)+i]=lln[x*(lastleaf+1)+j];
                            llt[x*(lastleaf+1)+i]=llt[x*(lastleaf+1)+j];
                        }
                    }
                    // right
                    if (rln[x*(lastleaf+1)+j]>0){
                        if (rln[x*(lastleaf+1)+i]!=0){
                            rlc[x*(lastleaf+1)+rlt[x*(lastleaf+1)+i]]=rlc[x*(lastleaf+1)+j];
                            rln[x*(lastleaf+1)+i]+=rln[x*(lastleaf+1)+j];
                            rlt[x*(lastleaf+1)+i]=rlt[x*(lastleaf+1)+j];
                        }else{
                            rlc[x*(lastleaf+1)+i]=rlc[x*(lastleaf+1)+j];
                            rln[x*(lastleaf+1)+i]=rln[x*(lastleaf+1)+j];
                            rlt[x*(lastleaf+1)+i]=rlt[x*(lastleaf+1)+j];
                        }
                    }
                }

                j=sst->sibling(j);
            }


        }

    }

//    // debug
//    for (int i=0; i<(int)pt.size(); i++){
//        if (len[i]>=10){
//            cout<<"("<<pt[i]<<","<<pq[i]<<","<<len[i]<<")"<<endl;
//        }
//    }

    // release memorey
    delete text;
    delete sst;
    delete llc;
    delete rlc;
    delete lln;
    delete rln;
    delete llt;
    delete rlt;
    delete depth;
}

/**
 * @brief NucleotideAlignmentMethod::exact_match_segment_chain
 * @param ems
 * @param si
 */
void NucleotideAlignmentMethod::exact_match_segment_chain(vector<X> &ems, vector<int> &si){
    if (ems.empty()) return;

    // comparison proxy
    struct Larger{
        bool operator()(const X& a, const X& b){
            return a.len>b.len;
        }
    }larger;

    // temporary array
    vector<X> B(ems);
    sort(B.begin(), B.end(), larger);

    // index of the largest exact match segment
    int i=B[0].i;

    // left-side and right-side
    vector<X> left, right;
    for (int t=0; t<(int)ems.size(); t++){
        if (ems[t].pt+ems[t].len-1<=B[0].pt && ems[t].pq+ems[t].len-1<=B[0].pq){
            left.push_back(ems[t]);
        }
        if (ems[t].pt>=B[0].pt+B[0].len-1 && ems[t].pq>=B[0].pq+B[0].len-1){
            right.push_back(ems[t]);
        }
    }

    // find on left side
    vector<int> lsi;
    exact_match_segment_chain(left, lsi);

    // find on right side
    vector<int> rsi;
    exact_match_segment_chain(right, rsi);

    // package
    for (int t=0; t<(int)lsi.size(); t++){
        si.push_back(lsi[t]);
    }
    si.push_back(i);
    for (int t=0; t<(int)rsi.size(); t++){
        si.push_back(rsi[t]);
    }
}

/**
 * @brief NucleotideAlignmentMethod::find_exact_match_segment_chain
 * @param align
 * @param t0
 * @param t1
 * @param q0
 * @param q1
 */
void NucleotideAlignmentMethod::find_exact_match_segment_chain(const NucleotideAlignment &align, int &t0, int &t1, int &q0, int &q1){
    vector<int> pt, pq, len;
    maximal_pairs(align, pt, pq, len);

    vector<X> ems;
    for (int i=0; i<(int)pt.size(); i++){
        X x(i, pt[i], pq[i], len[i]);
        ems.push_back(x);
    }

    vector<int> si;
    exact_match_segment_chain(ems, si);

    t0=ems[si[0]].pt;
    t1=ems[si[si.size()-1]].pt+ems[si[si.size()-1]].len-1;

    q0=ems[si[0]].pq;
    q1=ems[si[si.size()-1]].pq+ems[si[si.size()-1]].len-1;
}
