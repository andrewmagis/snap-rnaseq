#pragma once

//System includes
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <pthread.h>

#include "IntervalTree.h"
#include "Genome.h"
#include "FASTA.h"

#define GTF_MAX_READ_SIZE 4096

//Definitions for the possibilities of paired-end alignments          
const int FIRST_NOT_ALIGNED         = 0; 
const int SECOND_NOT_ALIGNED        = 1;
const int NOT_REVERSE_COMPLIMENTED  = 2;
const int ALIGNED_SAME_GENE         = 3;
const int ALIGNED_DIFF_GENE         = 4;
const int ALIGNED_SAME_CHR          = 5;
const int ALIGNED_DIFF_CHR          = 6;
const int CIRCULAR                  = 7;

class GTFReader;

//Namespaces
using namespace std;

class ReadInterval {

    friend class GTFReader;
    friend class ReadIntervalPair;
    friend class ReadIntervalMap;
    
    public:
        ReadInterval() : start(0), end(0), is_spliced(false), consolidated(false) {};
        ReadInterval(string chr, unsigned start_, unsigned end_, string id, bool is_spliced_);
        ReadInterval(const ReadInterval &rhs);
        virtual ~ReadInterval();
        ReadInterval& operator=(const ReadInterval &rhs);
        bool operator<(const ReadInterval &rhs) const;
        
        std::string Chr() const { return chr; };
        unsigned Start() const { return start; };
        unsigned End() const { return end; };
        std::string GeneID() const;
        std::string GeneName() const;
        std::string GeneNameSpliced(unsigned intersection) const;
        unsigned Intersection(const ReadInterval *rhs, set<string> &difference) const;
        unsigned Difference(const ReadInterval *rhs, set<string> &difference) const;
        void Print() const;
        void PrintIDs() const;
        void PrintWithMate();
        
        void WriteGTF(ofstream &outfile, unsigned intersection) const;
        void GetGeneInfo(GTFReader *gtf);
    
    protected:
    
        void UpdateMatePointers(ReadInterval *new_pointer);
    
        string chr;
        unsigned start;
        unsigned end;
        set<string> ids;
        set<string> gene_ids;
        set<string> gene_names;
        bool is_spliced;
        bool consolidated;
        
        //Pointer to mate(s)
        set<ReadInterval*> mate;

};

class ReadIntervalPair {

    friend class GTFReader;
    friend class ReadIntervalMap;
    
    public:
        ReadIntervalPair(ReadInterval *interval1, ReadInterval *interval2);
        ReadIntervalPair(const ReadIntervalPair &rhs);
        virtual ~ReadIntervalPair();
        ReadIntervalPair& operator=(const ReadIntervalPair &rhs);
        bool operator<(const ReadIntervalPair &rhs) const;
        unsigned Intersection() const { return intersection.size(); };
        
        void Write(ofstream &outfile) const;
        void WriteGTF(ofstream &outfile) const;
        void Print() const;
        void PrintIDs() const;

    protected:
        
        ReadInterval *interval1;
        ReadInterval *interval2;
        set<string> intersection;
        
};

typedef std::pair<ReadIntervalPair, ReadIntervalPair> interval_pair;
typedef std::vector<interval_pair> spliced_mate;

class ReadIntervalMap {

    public:
        ReadIntervalMap();
        ReadIntervalMap(const ReadIntervalMap &rhs);
        virtual ~ReadIntervalMap();
        ReadIntervalMap& operator=(const ReadIntervalMap &rhs);
        
        void AddInterval(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, bool is_spliced);
        void Consolidate(GTFReader *gtf, unsigned buffer);
        void Intersect(const ReadIntervalMap &map);
        void Print() const;
        void PrintSplicedMatePairs();
        
        void WriteSplicedMatePairs(ofstream &logfile, ofstream &readfile);
        void WriteGTF(ofstream &outfile);
        
    protected:

        unsigned ConsolidateReadIntervals(unsigned buffer);
        pthread_mutex_t mutex;
  
        //Vector of all paired-end reads that are not between genes
        std::vector<Interval<ReadInterval*> > read_intervals;
        IntervalTree<ReadInterval*> read_tree;  
        std::vector<ReadIntervalPair> pairs;  
        
        //Final results for a given map that contains overlaps between paired and spliced reads
        spliced_mate spliced_mate_pairs;

};

class GTFFeature {

    friend class GTFReader;
    friend class GTFTranscript;
    
    public:
    
        GTFFeature(string line);
        GTFFeature(const GTFFeature &rhs);
        virtual ~GTFFeature();
        GTFFeature& operator=(const GTFFeature &rhs);
        bool operator<(const GTFFeature &rhs) const;
        
        unsigned Length() const;
        bool GetAttribute(string key, string &value) const;
        void Print() const;
        
        
    protected:
    
        string chr;
        string source;
        string feature;
        unsigned start;
        unsigned end;
        string score;
        char strand;
        char frame;
        string gene_id;
        string transcript_id;
        string gene_name;
        std::map<std::string, std::string> attributes;

};

typedef std::vector<GTFFeature> feature_list;
typedef std::pair<unsigned, unsigned> junction;

class GTFTranscript {

    friend class GTFReader;
    
    public:
        
        GTFTranscript(string chr, string gene_id, string transcript_id, unsigned start, unsigned end);
        GTFTranscript(const GTFTranscript &rhs);
        virtual ~GTFTranscript();
        GTFTranscript& operator=(const GTFTranscript &rhs);
        
        string Chr() const { return chr; };
        string TranscriptID() const { return transcript_id; };
        string GeneID() const { return gene_id; };
        unsigned GenomicPosition(unsigned transcript_pos, unsigned span) const;
        void Junctions(unsigned start, unsigned span, std::vector<junction> &junctions) const;
        void WriteFASTA(const Genome *genome, std::ofstream &outfile) const;
        
    protected:
    
        void UpdateBoundaries(unsigned start, unsigned end);
        void Process();
    
        unsigned start;
        unsigned end;
        string chr;
        string gene_id;
        string transcript_id;
        feature_list features;
        feature_list exons;

};

class GTFGene {

    friend class GTFReader;
    
    public:
    
        GTFGene(string chr, string gene_id, unsigned start, unsigned end, string gene_name_);
        GTFGene(const GTFGene &rhs);
        virtual ~GTFGene();
        GTFGene & operator=(const GTFGene &rhs);
        bool operator<(const GTFGene &rhs) const;
        
        string Chr() const { return chr; };
        string GeneID() const { return gene_id; };
        string GeneName() const { return gene_name; };
        bool CheckBoundary(string query_chr, unsigned query_pos, unsigned buffer=1000) const;
        void Print() const;
        
    protected:
    
        void UpdateBoundaries(unsigned start, unsigned end);
    
        string chr;
        string gene_id;
        string gene_name;
        unsigned start;
        unsigned end;
                
};        

typedef std::map<std::string, GTFTranscript> transcript_map;
typedef std::map<std::string, GTFGene> gene_map;

class GTFReader {

    public:

        GTFReader();
        virtual ~GTFReader();
        
        int Load(string filename);
        const GTFTranscript& GetTranscript(string transcript_id) const;
        const GTFGene& GetGene(string gene_id) const;
        unsigned Size() const { return transcripts.size(); };
        void IntervalGenes(std::string chr, unsigned start, unsigned stop, std::vector<GTFGene> &results);
        
        //Functions for building the transcriptome file
        void BuildTranscriptome(const Genome *genome);
        
        void TransGenePair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id);
        void CisChromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id);
        void TransChromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id);
 
        void TransGeneSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id);
        void CisChromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id);
        void TransChromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id);
               
        void PrintGeneAssociations();
               
        void Test();
    
    protected:

        int Parse(string line);
        unsigned ConsolidateReadIntervals(unsigned buffer);
        
        string filename;
        transcript_map transcripts;
        gene_map genes;
        std::vector<Interval<GTFGene> > intervals;
        IntervalTree<GTFGene> tree;
        
        //Objects for read intervals
        ReadIntervalMap trans_gene_pairs;
        ReadIntervalMap cis_chromosomal_pairs;
        ReadIntervalMap trans_chromosomal_pairs;
        ReadIntervalMap trans_gene_splices;
        ReadIntervalMap cis_chromosomal_splices;
        ReadIntervalMap trans_chromosomal_splices;

};

inline std::string ToString(const unsigned& arg)
{
  char buffer[65];
  sprintf(buffer, "%u", arg) ;
  return std::string(buffer);
}

