/*++

Module Name:

    AlignmentFilter.h

Abstract:

    Filters transcriptome and genome alignments, allows for conversion from transcriptome
    to genome coordinates. Also handles fusion gene searching

Authors:

    Andrew Magis, June, 2013

Environment:

Revision History:

--*/

/*

#pragma once

//System includes
#include <fstream>
#include <string>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

#include "PairedEndAligner.h"
#include "SpecialAligner.h"
#include "Genome.h"
#include "GTFReader.h"

using namespace std;

class Alignment {

    friend class AlignmentFilter;
    friend class AlignmentPair;

    public:
    
        Alignment(unsigned location_, bool isRC_, int score_, string rname_, unsigned pos_, unsigned pos_end_, unsigned pos_original_, string transcript_id_, string gene_id_, bool isTranscriptome_);
        virtual ~Alignment() {}
        Alignment(const Alignment &rhs);
        Alignment& operator=(const Alignment&rhs);
        void Print();
        
    protected:
    
        unsigned location;
        bool isRC;
        unsigned score;
        string rname;
        unsigned pos;
        unsigned pos_end;
        unsigned pos_original;
        bool isTranscriptome;
        string transcript_id;
        string gene_id;
        string hashkey;
        
};

class AlignmentPair {

    friend class AlignmentFilter;

    public: 
    
        AlignmentPair(Alignment *align1_, Alignment *align2_, char flag_, bool is_unannotated_, bool is_backspliced_);
        virtual ~AlignmentPair() {}
        AlignmentPair(const AlignmentPair &rhs);
        AlignmentPair& operator=(const AlignmentPair &rhs);
        bool operator<(const AlignmentPair &rhs) const;
        void Print();
        
    protected:
    
        Alignment *align1;
        Alignment *align2;
        char flag;
        int distance;
        unsigned score;
        bool is_unannotated;
        bool is_backspliced;

};

typedef std::map<std::string, Alignment> alignment_map;

class AlignmentFilter {
  
    public:
    
        //Constructors/destructor
        AlignmentFilter(Read *read0_, Read *read1_, const Genome* genome_, const Genome* transcriptome_, GTFReader* gtf_, unsigned minSpacing_, unsigned maxSpacing_, unsigned confDiff_, unsigned maxDist_, unsigned seedLen_, SpecialAligner *specialAligner);
        virtual ~AlignmentFilter();  
        
        //Functions
        int AddAlignment(unsigned location, bool isRC, int score, bool isTranscriptome, bool isMate0); 
        int Filter(PairedAlignmentResult* result);    
        
    protected:
    
        int HashAlignment(Alignment& alignment, alignment_map& hashtable);
        void ConvertToGenomic();
        void ProcessPairs(PairedAlignmentResult* result, std::vector<AlignmentPair> &pairs);
        void CheckNoRC(PairedAlignmentResult* result, std::vector<AlignmentPair> &pairs);
        void FindPartialMatches(PairedAlignmentResult* result, AlignmentPair &pair);
        
        //Called on all unaligned reads to look for unknown splice junctions
        void UnalignedRead(Read *read, unsigned minDiff);
        bool ProcessSplices(std::vector<AlignmentPair> &pairs, unsigned minDiff);
              
        //Temp printing
        void PrintMaps(seed_map &map, seed_map &mapRC);
        

        Read *read0;
        Read *read1;
        const Genome* genome;
        const Genome* transcriptome;
        GTFReader* gtf;
        alignment_map mate0;
        alignment_map mate1;        
        unsigned minSpacing;
        unsigned maxSpacing;
        unsigned confDiff;
        unsigned maxDist;
        unsigned seedLen;
        SpecialAligner *specialAligner;

};

*/

