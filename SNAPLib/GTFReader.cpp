//Source includes
#include "GTFReader.h"

//System includes
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sstream>

#include "Compat.h"

bool IntervalNodeSort(const Interval<ReadInterval*> &interval0, const Interval<ReadInterval*> &interval1) { 
    return interval0.value->Start() < interval1.value->Start(); 
}

bool SplicedMateSort(const interval_pair &pair0, const interval_pair &pair1) {
    return ((pair0.first.Intersection() + pair0.second.Intersection()) >
            (pair1.first.Intersection() + pair1.second.Intersection()));

}

ReadInterval::ReadInterval(string chr_, unsigned start_, unsigned end_, string id, bool is_spliced_) 
    : chr(chr_), start(start_), end(end_), is_spliced(is_spliced_), consolidated(false)
{
    ids.insert(id);
} 
    
ReadInterval::ReadInterval(const ReadInterval &rhs) 
    : chr(rhs.chr), start(rhs.start), end(rhs.end), ids(rhs.ids), is_spliced(rhs.is_spliced), consolidated(rhs.consolidated), mate(rhs.mate)
{}

ReadInterval::~ReadInterval() 
{}

ReadInterval& ReadInterval::operator=(const ReadInterval &rhs) {
    if (this != &rhs) {
        chr = rhs.chr;
        start = rhs.start;
        end = rhs.end;
        ids = rhs.ids;
        is_spliced = rhs.is_spliced;
        consolidated = rhs.consolidated;
        mate = rhs.mate;
    }
    return *this;    
}

void ReadInterval::UpdateMatePointers(ReadInterval *new_pointer) {
    
    //Go through each mate
    for (set<ReadInterval*>::iterator it = mate.begin(); it != mate.end(); ++it) {
    
        //Delete the keys that point to this object
        (*it)->mate.erase(this); 
        
        //Replace it with the new pointer
        (*it)->mate.insert(new_pointer);
    
    }
}

unsigned ReadInterval::Intersection(const ReadInterval *rhs, set<string> &intersection) const {

    //Return the size of the intersection between the two sets of IDs
    for (set<string>::iterator it = ids.begin(); it != ids.end(); ++it) {
              
        set<string>::iterator pos = rhs->ids.find(*it);
        if (pos != rhs->ids.end()) {
            intersection.insert(*it);
        }        
    }
    return intersection.size();
}

unsigned ReadInterval::Difference(const ReadInterval *rhs, set<string> &difference) const {

    //Return the size of the intersection between the two sets of IDs
    for (set<string>::iterator it = ids.begin(); it != ids.end(); ++it) {
        
        set<string>::iterator pos = rhs->ids.find(*it);
        if (pos == rhs->ids.end()) {
            difference.insert(*it);
        }        
    }
    return difference.size();
}

bool ReadInterval::operator<(const ReadInterval& rhs) const {
    return (start <= rhs.start);
}

std::string ReadInterval::GeneID() const {

    std::set<string>::iterator it = gene_ids.begin();
    if (it == gene_ids.end()) {
        return "NoGene";
    }
    string gene_id = (*it);
    for (++it; it != gene_ids.end(); ++it) {
        gene_id += ',' + (*it);
    }
    return gene_id;
}

std::string ReadInterval::GeneName() const {

    std::set<string>::iterator it = gene_names.begin();
    if (it == gene_names.end()) {
        return GeneID();
    }
    string gene_name = (*it);
    for (++it; it != gene_names.end(); ++it) {
        gene_name += ',' + (*it);
    }
    return gene_name;
}

std::string ReadInterval::GeneNameSpliced(unsigned intersection) const {

    //Get gene_name, if present or gene_id if not
    string gene_name = GeneName();
    if (is_spliced) {
        gene_name += ",S";
    } else {
        gene_name += ",P";
    }
    return gene_name + "," + ToString(intersection);
}

void ReadInterval::GetGeneInfo(GTFReader *gtf) {

    //Get gene info for this interval
    std::vector<GTFGene> results;
    gtf->IntervalGenes(chr, start, end, results);
    
    for (std::vector<GTFGene>::iterator it = results.begin(); it != results.end(); ++it) {
        gene_ids.insert(it->GeneID());
        if (it->GeneName().size() > 0) {
            gene_names.insert(it->GeneName());
        }
    }

}

void ReadInterval::Write(ofstream &outfile, unsigned intersection) const {
    outfile << chr << ':' << start << '-' << end << '\t' << GeneID() << '\t' << GeneNameSpliced(intersection) << endl;
}

void ReadInterval::Print() const {
    cout << chr << ':' << start << '-' << end << '\t' << GeneID() << endl;
}

void ReadInterval::WriteGTF(ofstream &outfile, unsigned intersection) const {
    outfile << chr << '\t' << "snap-rna" << '\t' << "interval" << '\t' << start << '\t' << end << '\t' << '.' << '\t' << '.' << '\t' << '.' << '\t';
    outfile << "gene_id \"" << GeneID() << "\"; transcript_id \"" << GeneNameSpliced(intersection) << "\"; gene_name \"" << GeneName() << "\";" << endl;
}

ReadIntervalPair::ReadIntervalPair(ReadInterval *interval1_, ReadInterval *interval2_) 
    : interval1(interval1_), interval2(interval2_)
{
    interval1->Intersection(interval2, intersection);
} 
    
ReadIntervalPair::ReadIntervalPair(const ReadIntervalPair &rhs) 
    : interval1(rhs.interval1), interval2(rhs.interval2), intersection(rhs.intersection)
{}

ReadIntervalPair::~ReadIntervalPair() 
{}

ReadIntervalPair& ReadIntervalPair::operator=(const ReadIntervalPair &rhs) {
    if (this != &rhs) {
        interval1 = rhs.interval1;
        interval2 = rhs.interval2;
        intersection = rhs.intersection;
    }
    return *this;    
}

bool ReadIntervalPair::operator<(const ReadIntervalPair& rhs) const {
    return intersection.size() > rhs.intersection.size();
}

void ReadIntervalPair::Write(ofstream &outfile) const {
    interval1->Write(outfile, intersection.size());
    interval2->Write(outfile, intersection.size());
}

void ReadIntervalPair::WriteGTF(ofstream &outfile) const {
    interval1->WriteGTF(outfile, intersection.size());
    interval2->WriteGTF(outfile, intersection.size());
}

ReadIntervalMap::ReadIntervalMap() {   
    pthread_mutex_init(&mutex, NULL);
}
   
ReadIntervalMap::ReadIntervalMap(const ReadIntervalMap &rhs) {
    printf("ReadIntervalMap: copy constructor not implemented\n");
    exit(1);
}

ReadIntervalMap::~ReadIntervalMap() {
    pthread_mutex_destroy(&mutex);

    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
        delete it->value;
    }
    read_intervals.clear();
}

void ReadIntervalMap::Clear() {
    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
        delete it->value;
    }
    read_intervals.clear();
}

ReadIntervalMap& ReadIntervalMap::operator=(const ReadIntervalMap &rhs) {
    printf("ReadIntervalMap: assignment operator not implemented\n");
    exit(1);
}

void ReadIntervalMap::AddInterval(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id, bool is_spliced) {

    //Get the mutex
    pthread_mutex_lock(&mutex);
    
    ReadInterval *mate0 = new ReadInterval(chr0, start0, end0, id, is_spliced);
    ReadInterval *mate1 = new ReadInterval(chr1, start1, end1, id, is_spliced);

    mate0->mate.insert(mate1);
    mate1->mate.insert(mate0);
    
    //Add this to the vector of links
    read_intervals.push_back(Interval<ReadInterval*>(start0, end0, mate0));
    read_intervals.push_back(Interval<ReadInterval*>(start1, end1, mate1));
    
    //Unlock it
    pthread_mutex_unlock(&mutex);
    
}

void ReadIntervalMap::Consolidate(GTFReader *gtf, unsigned buffer=100) {

    printf("Building interval tree of read pairs\n");
    
    //Get size of current interval set
    unsigned initial_size = read_intervals.size();
       
    //Repeatedly consolidate existing reads until it no longer changes
    do {
        initial_size = read_intervals.size();
        ConsolidateReadIntervals(buffer);
        printf("Initial Size: %u Current Size: %u\n", initial_size, read_intervals.size());
    } while (initial_size > read_intervals.size());

    
    //Now we build all ReadIntervalPairs
    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
        
        //Get interval info 
        it->value->GetGeneInfo(gtf);
        
        //For each mate in this result
        for (std::set<ReadInterval*>::iterator mate_it = it->value->mate.begin(); mate_it != it->value->mate.end(); ++mate_it) {    
            pairs.push_back(ReadIntervalPair(it->value, *mate_it));      
        }
    }
    
    //Sort the pairs by the number of mate pairs they have in common
    sort(pairs.begin(), pairs.end());    
    
}

unsigned ReadIntervalMap::ConsolidateReadIntervals(unsigned buffer) {

    //Sort the read intervals, which avoids weird issue in IntervalTree
    sort(read_intervals.begin(), read_intervals.end(), IntervalNodeSort);

    //Build the interval tree from the current set of intervals
    read_tree = IntervalTree<ReadInterval*>(read_intervals);
    
    //Create a new, temporary set of intervals
    std::vector<Interval<ReadInterval*> > temp_intervals;
    
    //Loop over all the read_intervals that we currently have
    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {

        //Check to see if start is 0, which means it has already been consolidated
        if (it->value->consolidated) {
            continue;
        }

        //Query the interval tree for this interval to find any overlaps
        std::vector<Interval<ReadInterval*> > temp_results;
        std::vector<ReadInterval*> chr_results;
        read_tree.findOverlapping(it->value->start-buffer, it->value->end+buffer, temp_results);
 
        //Filter by chromosome
        for (std::vector<Interval<ReadInterval*> >::iterator it2 = temp_results.begin(); it2 != temp_results.end(); ++it2) {
        
            //Make sure this has not already been consolidated
            if (it2->value->consolidated) {
                continue;
            }
        
            //Make sure the chromosomes match
            if (it2->value->chr.compare(it->value->chr) == 0) {
                chr_results.push_back(it2->value);
            }   
        }   

        //I have already been consolidated
        if (chr_results.size() == 0) {
            printf("Warning, cannot find myself\n");
            continue;
        } 
        
        //In this case this read has no overlaps
        if (chr_results.size() == 1) {
            //Add this read back into the result set 
            temp_intervals.push_back(Interval<ReadInterval*>(chr_results[0]->start, chr_results[0]->end, chr_results[0]));
            continue;
        }
              
        //Now we have all the overlapping reads of this read +/- buffer in the vector chr_results
        ReadInterval *new_interval = new ReadInterval;
        for (std::vector<ReadInterval*>::iterator it3 = chr_results.begin(); it3 != chr_results.end(); ++it3) {
             
            //Build a new interval with the existing read IDs    
            if (new_interval->start == 0) {
                new_interval->chr = (*it3)->chr;
                new_interval->start = (*it3)->start;
                new_interval->end = (*it3)->end;
                new_interval->ids = (*it3)->ids;
                new_interval->is_spliced = (*it3)->is_spliced;

            //Otherwise extend the boundary of this interval
            } else {
            
                new_interval->start = min(new_interval->start, (*it3)->start);
                new_interval->end = max(new_interval->end, (*it3)->end);
                std::copy((*it3)->ids.begin(), (*it3)->ids.end(), std::inserter(new_interval->ids, new_interval->ids.end()));
       
            }            
            
            //Update the boundaries of any intervals that pointed to this interval
            (*it3)->UpdateMatePointers(new_interval);
            
            //Update the mates
            for (set<ReadInterval*>::iterator it4 = (*it3)->mate.begin(); it4 != (*it3)->mate.end(); ++it4) {
                new_interval->mate.insert(*it4);
            }
            
            //Indicate this interval has been used
            (*it3)->consolidated = true;
            
        }
                
        //Now we have built new interval. Insert it in the new interval vector
        temp_intervals.push_back(Interval<ReadInterval*>(new_interval->start, new_interval->end, new_interval));
    
    } 
    
    /*
    //Delete all the original intervals that were zeroed out
    for (std::vector<Interval<ReadInterval*> >::iterator it = read_intervals.begin(); it != read_intervals.end(); ++it) {
        if (it->value->consolidated) {
            delete it->value;
        }
    }
    */
    
    //Finally, assign set of new intervals to the old set, and return the size
    read_intervals = temp_intervals;
    return read_intervals.size();

}

void ReadIntervalMap::Intersect(const ReadIntervalMap &rhs, unsigned buffer=100) {

    std::vector<Interval<ReadInterval*> > left_intervals;
    std::vector<Interval<ReadInterval*> > right_intervals;

    //Build the new interval tree
    for (std::vector<ReadIntervalPair>::iterator it = pairs.begin(); it != pairs.end(); ++it) {

        left_intervals.push_back(Interval<ReadInterval*>(it->interval1->start, it->interval1->end, it->interval1));
        right_intervals.push_back(Interval<ReadInterval*>(it->interval2->start, it->interval2->end, it->interval2));

    }
    
    //Sort these intervals by start position
    sort(left_intervals.begin(), left_intervals.end(), IntervalNodeSort);  
    sort(right_intervals.begin(), right_intervals.end(), IntervalNodeSort);
    
    IntervalTree<ReadInterval*> left_tree(left_intervals);
    IntervalTree<ReadInterval*> right_tree(right_intervals);
    
    //Loop over all intervals in the other map
    for (std::vector<ReadIntervalPair>::const_iterator it = rhs.pairs.begin(); it != rhs.pairs.end(); ++it) {
    
        std::vector<Interval<ReadInterval*> > left_temp_results;
        std::vector<ReadInterval*> left_chr_results;
        read_tree.findOverlapping(it->interval1->start-buffer, it->interval1->end+buffer, left_temp_results);
  
        //Filter by chromosome
        for (std::vector<Interval<ReadInterval*> >::iterator it2 = left_temp_results.begin(); it2 != left_temp_results.end(); ++it2) {
    
            //Make sure the chromosomes match
            if (it2->value->chr.compare(it->interval1->chr) == 0) {
                left_chr_results.push_back(it2->value);
            }   
        } 
        
        std::vector<Interval<ReadInterval*> > right_temp_results;
        std::vector<ReadInterval*> right_chr_results;
        read_tree.findOverlapping(it->interval2->start-buffer, it->interval2->end+buffer, right_temp_results);
         
        //Filter by chromosome
        for (std::vector<Interval<ReadInterval*> >::iterator it2 = right_temp_results.begin(); it2 != right_temp_results.end(); ++it2) {
    
            //Make sure the chromosomes match
            if (it2->value->chr.compare(it->interval2->chr) == 0) {
                right_chr_results.push_back(it2->value);
            }   
        }     
        
        //Now all we have to do is verify that the results are linked to each other
        for (std::vector<ReadInterval*>::iterator left = left_chr_results.begin(); left != left_chr_results.end(); ++left) {
            for (std::vector<ReadInterval*>::iterator right = right_chr_results.begin(); right != right_chr_results.end(); ++right) {
        
                set<ReadInterval*>::iterator pos = (*left)->mate.find(*right);
        
                //If the mate of one contains the other
                if ((pos != (*left)->mate.end())) {

                    //Create a new ReadIntervalPair with both of these overlapping sets
                    spliced_mate_pairs.push_back( std::make_pair(ReadIntervalPair(it->interval1, it->interval2), ReadIntervalPair(*left, *right)) );

                }
            }
        }
    }
}

void ReadIntervalMap::WriteGTF(ofstream &outfile) {

    //Sort the pairs
    sort(spliced_mate_pairs.begin(), spliced_mate_pairs.end(), SplicedMateSort);

    //Write out each one to the output file
    for (spliced_mate::const_iterator it = spliced_mate_pairs.begin(); it != spliced_mate_pairs.end(); ++it) {
    
        it->first.WriteGTF(outfile);
        it->second.WriteGTF(outfile);
       
    }
}

void ReadIntervalMap::WriteSplicedMatePairs(ofstream &logfile, ofstream &readfile) {

    //Sort the pairs
    sort(spliced_mate_pairs.begin(), spliced_mate_pairs.end(), SplicedMateSort);

    //Write out each one to the output file
    for (spliced_mate::const_iterator it = spliced_mate_pairs.begin(); it != spliced_mate_pairs.end(); ++it) {
     
        it->first.Write(logfile);
        it->second.Write(logfile);
        logfile << endl;
        
    }
}

GTFFeature::GTFFeature(string line) 
    : count(0) 
{
                
    char* line_c = (char*)line.c_str(); 
    char *pch;
    pch = strtok(line_c,"'\t'"); chr = pch;
    pch = strtok(NULL,"'\t'"); source = pch;
    pch = strtok(NULL,"'\t'"); feature = pch;
    pch = strtok(NULL,"'\t'"); start = atoi(pch);
    pch = strtok(NULL,"'\t'"); end = atoi(pch);
    pch = strtok(NULL,"'\t'"); score = pch;
    pch = strtok(NULL,"'\t'"); strand = *pch;
    pch = strtok(NULL,"'\t'"); frame = *pch;
        
    while (pch != NULL) {
    
        char *temp_key, *temp_value;
        pch = strtok(NULL, " ="); temp_key = pch;
        pch = strtok(NULL, ";"); temp_value = pch;
        
        if (temp_key != NULL) {
            
            string key = temp_key;
            string value = temp_value;
            
            //Remove any quotes or spaces from each string
            value.erase(remove(value.begin(), value.end(), '\"' ), value.end());
            //printf("value: %s key: %s\n", value.c_str(), key.c_str());
            attributes.insert(std::map<string,string>::value_type(key, value));
        }
    }
    
    string value;
    
    //If gene_id is present, use this as gene_id
    if (GetAttribute("gene_id", value)) {
        gene_id = value;
        
    //Otherwise, use Parent
    } else if (GetAttribute("Parent", value)) {
        gene_id = value;
    
    } else {
        //printf("Cannot find gene identifiers 'gene_id' (GTF) or 'Parent' (GFF3) in annotation file\n");
        //exit(1);
        gene_id = "gene_id";
    } 
    
    //If transcript_id exists
    if (GetAttribute("transcript_id", value)) {
        transcript_id = value;
    } else {
        transcript_id = gene_id;
    }
    
    //if gene_name exists
    if (GetAttribute("gene_name", value)) {
        gene_name = value;
    } else if (GetAttribute("Name", value)) {
        gene_name = value;
    } else {
        gene_name = "";
    }
    
    
}

bool GTFFeature::GetAttribute(string key, string &value) const {

    //At the end we must define the gene_id and transcript_id
    std::map<string,string>::const_iterator pos = attributes.find(key);
    if (pos == attributes.end()) {
        return false;
    }
    
    //Return the attribute
    value = pos->second;
    return true;
}

unsigned GTFFeature::Length() const {
    return (end-start)+1;
}

void GTFFeature::Print() const {
    printf("%s\t%s\t%s\t%u\t%u\t%s\t%c\t%c\t%s\t%s\n", chr.c_str(), source.c_str(), feature.c_str(), start, end, score.c_str(), strand, frame, gene_id.c_str(), transcript_id.c_str());
}

GTFFeature::GTFFeature(const GTFFeature& rhs) 
    :   chr(rhs.chr), source(rhs.source), feature(rhs.feature), start(rhs.start), end(rhs.end), score(rhs.score),
        strand(rhs.strand), frame(rhs.frame), gene_id(rhs.gene_id), transcript_id(rhs.transcript_id), count(rhs.count)
{}
  
GTFFeature::~GTFFeature() {}

GTFFeature& GTFFeature::operator=(const GTFFeature& rhs) {

    if (this != &rhs) {
        chr = rhs.chr;
        source = rhs.source;
        feature = rhs.feature;
        start = rhs.start;
        end = rhs.end;
        score = rhs.score;
        strand = rhs.strand;
        frame = rhs.frame;
        gene_id = rhs.gene_id;
        transcript_id = rhs.transcript_id;
        count = rhs.count;
    }
    return *this;
}

bool GTFFeature::operator<(const GTFFeature& rhs) const {
    return start < rhs.start;
}

GTFGene::GTFGene(string _chr, string _gene_id, unsigned _start, unsigned _end, string gene_name_) 
    : chr(_chr), gene_id(_gene_id), start(_start), end(_end), gene_name(gene_name_)
{}

GTFGene::GTFGene(const GTFGene& rhs) 
    : chr(rhs.chr), gene_id(rhs.gene_id), start(rhs.start), end(rhs.end), gene_name(rhs.gene_name)
{}
  
GTFGene::~GTFGene() {}

GTFGene& GTFGene::operator=(const GTFGene& rhs) {

    if (this != &rhs) {
        chr = rhs.chr;
        gene_id = rhs.gene_id;
        start = rhs.start;
        end = rhs.end;
        gene_name = rhs.gene_name;
    }
    return *this;
}

bool GTFGene::operator<(const GTFGene &rhs) const {
    printf("< operator not implemented\n");
    exit(0);
}

void GTFGene::UpdateBoundaries(unsigned new_start, unsigned new_end) {
    start = std::min(new_start, start);
    end = std::max(new_end, end);
}   

bool GTFGene::CheckBoundary(string query_chr, unsigned query_pos, unsigned buffer) const {

    //Compare chromosomes
    if (chr.compare(query_chr) != 0) {
        return false;
    }
  
    //Compare position within buffered gene coordinates
    if ((query_pos >= std::max(start-buffer+1, (unsigned)1)) && (query_pos <= end+buffer)) {
        return true;
    }
    return false;
}

void GTFGene::Print() const {
    printf("%s\t%u\t%u\t%s\n", chr.c_str(), start, end, gene_id.c_str());
}

GTFTranscript::GTFTranscript(string _chr, string _gene_id, string _transcript_id, unsigned _start, unsigned _end) 
    : chr(_chr), gene_id(_gene_id), transcript_id(_transcript_id), start(_start), end(_end)
{}

GTFTranscript::GTFTranscript(const GTFTranscript& rhs) 
    : chr(rhs.chr), gene_id(rhs.gene_id), transcript_id(rhs.transcript_id), features(rhs.features), start(rhs.start), end(rhs.end)
{}
  
GTFTranscript::~GTFTranscript() {}

GTFTranscript& GTFTranscript::operator=(const GTFTranscript& rhs) {

    if (this != &rhs) {
        chr = rhs.chr;
        gene_id = rhs.gene_id;
        transcript_id = rhs.transcript_id;
        features = rhs.features;   
        start = rhs.start;
        end = rhs.end;
    }
    return *this;
}

void GTFTranscript::Process() {

    //Copy all exons into separate vector
    for (feature_list::iterator it = features.begin(); it != features.end(); ++it) {
        if (it->feature.compare("exon") == 0) {
            exons.push_back(*it);
        }
    }

    //Sort the exons by start and end
    sort(exons.begin(), exons.end());
    
}

void GTFTranscript::UpdateBoundaries(unsigned new_start, unsigned new_end) {
    start = std::min(new_start, start);
    end = std::max(new_end, end);
}   

unsigned GTFTranscript::GenomicPosition(unsigned transcript_pos, unsigned span) const {
              
    //This assumes 1-offset transcript pos, and returns 1-offset genomic position
    //Converts transcript coordinates into genomic coordinates
    for (feature_list::const_iterator it = exons.begin(); it != exons.end(); ++it) {
                
        //If transcript_pos is less than or equal to this
        if (transcript_pos > it->Length()) {
            transcript_pos -= it->Length();
        } else {
        
            //Get current position for start of alignment
            unsigned genome_pos = it->start + transcript_pos - 1;
                
            //Check to see if read exceeds length of transcript.  If so, return 0
            //This is possible due to the way SNAP uses consecutive chromosomes.
            if (genome_pos+span > end) {
                //printf("genome_pos: %u span: %u end: %u\n", genome_pos, span, end);
                //printf("Warning, transcript_pos exceeds transcript length\n");
                return 0;
                
            }
            return genome_pos;
        }
    }
    //printf("Warning, transcript_pos exceeds transcript length\n");
    return 0;
}

void GTFTranscript::Junctions(unsigned transcript_pos, unsigned span, std::vector<junction> &junctions) const {
    
    //printf("Transcript: %s %u %u\n", transcript_id.c_str(), transcript_pos, span);
    
    //Go through each feature of this transcript until we find the start
    unsigned current_pos = 0;
    feature_list::const_iterator current = exons.begin();
    for (feature_list::const_iterator next = ++(exons.begin()); next != exons.end(); ++next) {

        //Get the end of this exon
        current_pos += current->Length();
        
        //printf("Transcript Pos: %u Current Pos: %u \n", transcript_pos,current_pos);
        
        //If transcript_pos is less than or equal to current, the start
        //lies within this feature
        if (transcript_pos <= current_pos) {

            //printf("Inside: %u - %u + 1 >= %u\n", current_pos, transcript_pos, span);

            //Check to see if the span exceeds the current feature. If so, 
            //simply return the current list of junctions
            if (current_pos-transcript_pos+1 >= span) {
                return;
                
            //Otherwise, add the junction to the list of junctions
            } else {
            
                junctions.push_back(junction(current_pos+1, next->start-current->end-1));
                span -= (current_pos-transcript_pos+1);
                transcript_pos += (current_pos-transcript_pos+1);  
            
            }
        }
        current = next;         
    }

}

void GTFTranscript::WriteFASTA(const Genome *genome, std::ofstream &outfile) const {

    //Get the offset for this chromosome
    unsigned offset, amountExceeded;
    bool isValid = genome->getOffsetOfPiece(chr.c_str(), &offset);

    if (!isValid) {
        printf("Warning: chromosome %s from the annotation is not found in the genome file\n", chr.c_str());
        return;
    }
    
    string sequence;
    for (feature_list::const_iterator it = exons.begin(); it != exons.end(); ++it) {
    
        //Get the sequence of this feature from the genome
        const char* seq = genome->getSubstring(it->start+offset-1,  it->Length(), amountExceeded);
        
        if (seq == NULL) {
            printf("Warning: transcript %s exceeds boundaries of chromosome %s. Truncating\n", transcript_id.c_str(), chr.c_str());
            const char* seq = genome->getSubstring(it->start+offset-1, it->Length()-amountExceeded, amountExceeded);
            exit(1);
        }
        
        //Temp copy into string variable
        string temp(seq, it->Length());
        sequence += temp;
    
    }
    outfile << ">" + transcript_id << endl;
    outfile << sequence << endl;

}

//Constructor
GTFReader::GTFReader() {}

//Destructor
GTFReader::~GTFReader() {}

int GTFReader::Load(string _filename) {

    //Save the input filename
    filename = _filename;
    
    printf("Loading genome annotation from file... ");
    fflush(stdout);
    _int64 loadStart = timeInMillis();
    
    //Open input file
    std::ifstream infile(filename.c_str(), std::ios::in);
    if (!infile.is_open()) { 
        printf("Cannot open input file '%s'\n", filename.c_str());
	    exit(1);
    }
    
    string line;
    unsigned num_lines = 0;
    std::getline(infile, line, '\n');
    while (!infile.eof()){
  
        //Process this line
        Parse(line);
        num_lines++;
          
        std::getline(infile, line, '\n');
    }
    infile.close();
    
    //Sort each transcript
    for (transcript_map::iterator it = transcripts.begin(); it != transcripts.end(); ++it) {
        (*it).second.Process();
    }
    
    //Add all genes to interval tree
    for (gene_map::iterator it = genes.begin(); it != genes.end(); ++it) {
        gene_intervals.push_back(Interval<GTFGene*>(it->second.start, it->second.end, &it->second));
    }
    gene_tree = IntervalTree<GTFGene*>(gene_intervals);
    
    //Add all transcripts to interval tree
    for (transcript_map::iterator it = transcripts.begin(); it != transcripts.end(); ++it) {
        transcript_intervals.push_back(Interval<GTFTranscript*>(it->second.start, it->second.end, &it->second));
    }
    transcript_tree = IntervalTree<GTFTranscript*>(transcript_intervals);    
    
    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds. %u features, %u transcripts\n", loadTime / 1000, num_lines, transcripts.size());
    
}

int GTFReader::Parse(string line) {

    //Check to see if this is a comment
    if (line[0] == '#') {
        return 1;
    }

    // Create a new GTFFeature from this line
    GTFFeature feature(line);
    
    if (feature.feature.compare("exon") != 0) {
        return 1;
    }
    
    //We don't try to find features, we just add them to the vector and to the tree
    features.push_back(feature);
    
    //Try to find this transcript in the transcript_map
    transcript_map::iterator pos = transcripts.find(feature.transcript_id);
        
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    if ((pos == transcripts.end())) {
        GTFTranscript transcript(feature.chr, feature.gene_id, feature.transcript_id, feature.start, feature.end);
        transcript.features.push_back(feature);
        transcripts.insert(transcript_map::value_type(feature.transcript_id, transcript));
        
    //Otherwise, add this feature to the transcript
    } else {
        pos->second.features.push_back(feature);
        pos->second.UpdateBoundaries(feature.start, feature.end);
    }
    
    
    //Try to find this gene in the gene_map
    gene_map::iterator gpos = genes.find(feature.gene_id);
    
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    if ((gpos == genes.end())) {
        GTFGene gene(feature.chr, feature.gene_id, feature.start, feature.end, feature.gene_name);
        //gene.features.push_back(feature);
        genes.insert(gene_map::value_type(feature.gene_id, gene));
        
    //Otherwise, add this feature to the transcript
    } else {
        //Do nothing
        gpos->second.UpdateBoundaries(feature.start, feature.end);
        //gspos->second.features.push_back(feature);
    }    
    
    
        
    return 0;
    
}

const GTFTranscript& GTFReader::GetTranscript(string transcript_id) const {
 
    transcript_map::const_iterator pos = transcripts.find(transcript_id);
    if (pos == transcripts.end()) {
        //raise exception
        printf("No transcript %s\n", transcript_id.c_str());
        exit(1);
    }
    return pos->second;

}

const GTFGene& GTFReader::GetGene(string gene_id) const {
 
    gene_map::const_iterator pos = genes.find(gene_id);
    if (pos == genes.end()) {
        //raise exception
        printf("No gene %s\n", gene_id.c_str());
        exit(1);
    }
    return pos->second;

}

void GTFReader::IncrementReadCount(string transcript_id0, unsigned start0, unsigned end0, string transcript_id1, unsigned start1, unsigned end1) {

    //If a read is aligned to transcriptome, it can be spliced
    //If it is not aligned to transcriptome, it cannot be spliced
    
    //If the first read is aligned to transcriptome
    if (transcript_id0.size() != 0) {
    
        //Get the transcript associated with this id
        const GTFTranscript &transcript0 = GetTranscript(transcript_id0);
        
        //Get the junctions for this transcript
        std::vector<junction> junctions;
        transcript0.Junctions(start0, end0-start0+1, junctions);
        
        //For each junction, query the exon tree
            
    
    }

    //Get the transcript for this gene
    if (transcript_id1.size() != 0) {
    
        const GTFTranscript &transcript1 = GetTranscript(transcript_id1);
    
    
    }
    
    
    
}

void GTFReader::IntervalGenes(std::string chr, unsigned start, unsigned stop, std::vector<GTFGene> &results) {
    
    std::vector<Interval<GTFGene*> > temp_results;
    gene_tree.findOverlapping(start, stop, temp_results);

    //Filter by chromosome
    for (std::vector<Interval<GTFGene*> >::iterator it = temp_results.begin(); it != temp_results.end(); ++it) {
        if (it->value->chr.compare(chr) == 0) {
            results.push_back(*(it->value));
        }
    }
}

void GTFReader::IntrageneUnannotatedPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    intragene_unannotated_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, false);
}

void GTFReader::IntrageneUnannotatedSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    intragene_unannotated_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id, true);
}

void GTFReader::IntrageneCircularPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    intragene_circular_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, false);
}

void GTFReader::IntrageneCircularSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    intragene_circular_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id, true);
}

void GTFReader::IntrachromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    intrachromosomal_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, false);
}

void GTFReader::IntrachromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    intrachromosomal_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id, true);
}

void GTFReader::InterchromosomalPair(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    interchromosomal_pairs.AddInterval(chr0, start0, end0, chr1, start1, end1, id, false);
}

void GTFReader::InterchromosomalSplice(string chr0, unsigned start0, unsigned end0, string chr1, unsigned start1, unsigned end1, string id) {
    interchromosomal_splices.AddInterval(chr0, start0, end0, chr1, start1, end1, id, true);
}


void GTFReader::PrintGeneAssociations() {

    //Open output file
    ofstream interchromosomal_file, intrachromosomal_file, unannotated_file, circular_file;
    ofstream logfile, readfile;
    interchromosomal_file.open("interchromosomal_intervals.gtf");
    intrachromosomal_file.open("intrachromosomal_intervals.gtf");
    unannotated_file.open("unannotated_intervals.gtf");
    circular_file.open("circular_intervals.gtf");
	logfile.open("read_intervals.txt");
	readfile.open("read_ids.txt");

    interchromosomal_pairs.Consolidate(this, 100);
    interchromosomal_splices.Consolidate(this, 0);
    interchromosomal_pairs.Intersect(interchromosomal_splices, 100);
    logfile << "Inter-Chromosomal Intervals" << endl;
    interchromosomal_pairs.WriteGTF(interchromosomal_file);
    interchromosomal_pairs.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    interchromosomal_pairs.Clear();
    
    intrachromosomal_pairs.Consolidate(this, 100);
    intrachromosomal_splices.Consolidate(this, 0);
    intrachromosomal_pairs.Intersect(intrachromosomal_splices, 100);
    logfile << "Intra-Chromosomal Intervals" << endl;
    intrachromosomal_pairs.WriteGTF(intrachromosomal_file);
    intrachromosomal_pairs.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    intrachromosomal_pairs.Clear();
    
    intragene_unannotated_pairs.Consolidate(this, 100);
    intragene_unannotated_splices.Consolidate(this, 0);
    intragene_unannotated_pairs.Intersect(intragene_unannotated_splices, 100);
    logfile << "Intra-Gene Unannotated Intervals" << endl;
    intragene_unannotated_pairs.WriteGTF(unannotated_file);
    intragene_unannotated_pairs.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    intragene_unannotated_pairs.Clear();
    
    intragene_circular_pairs.Consolidate(this, 100);
    intragene_circular_splices.Consolidate(this, 0);
    intragene_circular_pairs.Intersect(intragene_circular_splices, 100); 
    logfile << "Intra-Gene Circular Intervals" << endl;
    intragene_circular_pairs.WriteGTF(circular_file);
    intragene_circular_pairs.WriteSplicedMatePairs(logfile, readfile);
    logfile << endl;
    intragene_circular_pairs.Clear();
    
    interchromosomal_file.close();
    intrachromosomal_file.close();
    unannotated_file.close();
    circular_file.close();
    logfile.close();
    readfile.close();
    
}

void GTFReader::Test() {

    /*
    GTFTranscript &transcript = GetTranscript("ENST00000489673");

    std::vector<junction> junctions = transcript.Junctions(1, 200);
    
    for (std::vector<junction>::iterator it = junctions.begin(); it != junctions.end(); ++it) {
        printf("[%d %d]\n", it->first, it->second);
    }
    */
    
    std::vector<GTFGene> results;
    IntervalGenes("22", 42290565, 42359490, results);
    printf("Size: %d\n", results.size());
    for (std::vector<GTFGene>::iterator it = results.begin(); it != results.end(); ++it) {
        it->Print();
    }


}

void GTFReader::BuildTranscriptome(const Genome *genome) {

    printf("Building transcriptome FASTA file... ");
    fflush(stdout);
    _int64 loadStart = timeInMillis();

    //Save the input filename
    string filename = "transcriptome.fa";
        
    //Open input file
    std::ofstream outfile(filename.c_str(), std::ios::out);
    if (!outfile.is_open()) { 
        printf("Cannot open output file '%s'\n", filename.c_str());
	    exit(1);
    }
    
    //Loop over all transcripts
    for (transcript_map::iterator it = transcripts.begin(); it != transcripts.end(); ++it) { 
        it->second.WriteFASTA(genome, outfile); 
    }
    
    //Close the output file
    outfile.close();

    _int64 loadTime = timeInMillis() - loadStart;
    printf("%llds.\n", loadTime / 1000);

}

/*
int main(int argc, char *argv[]) {

    const char *fasta_filename = NULL;
    const char *gtf_filename = NULL;

    //Process the command-line arguments
    int c;
    while ((c=getopt(argc, argv, "g:f:")) != EOF) {
    
        switch(c) {
            case 'g':
                gtf_filename = optarg;
                break;
            case 'f':
                fasta_filename = optarg;
            case ':':
                cerr << "Invalid option " << optarg << endl;
                exit(1);
                break; 	
	    }
    }
    
    struct timeval t1;
    gettimeofday(&t1, NULL);
    
    //Create the gtfReader
    GTFReader gtf;
    gtf.Load(gtf_filename);
    
    const Genome *genome = ReadFASTAGenome(fasta_filename);
    //gtf.Test();
    
    struct timeval t2;
    gettimeofday(&t2, NULL);
    
    double elapsedTime;
    elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
    elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;
    fprintf(stderr, "Time elapsed: %.3lf seconds\n", elapsedTime / 1000.0);
    return 0;
}

*/

/*
void GTFReader::AddCircRNA(string gene, bool isSpliced, const char* id, unsigned length) {

     //Try to find this key in the graph
    circ_graph::iterator pos = circ.find(gene);
    
    //Create a string for this read id
    string read_id(id, length);
        
    //If this sequence is not found, create a new vector to store this sequence (and others like it)
    std::pair<circ_graph::iterator,bool> ret;
    if ((pos == circ.end())) {
        ret = circ.insert(circ_graph::value_type(gene, gene));
        
        if (isSpliced) {
            ret.first->second.IncrementSplicedCount();
            ret.first->second.AddSplicedReadID(read_id);
        } else {
            ret.first->second.IncrementPairedCount();
            ret.first->second.AddPairedReadID(read_id);
        }
        
    //Otherwise, add this feature to the transcript
    } else {
        if (isSpliced) {
            pos->second.IncrementSplicedCount();
            pos->second.AddSplicedReadID(read_id);
        } else {
            pos->second.IncrementPairedCount();
            pos->second.AddPairedReadID(read_id);
        }
    }   

}*/

/*
GTFFeature::GTFFeature(string line) {

    char _chr[GTF_MAX_READ_SIZE];
    char _source[GTF_MAX_READ_SIZE];
    char _feature[GTF_MAX_READ_SIZE];
    unsigned _start;
    unsigned _end;
    char _score[GTF_MAX_READ_SIZE];
    char _strand;
    char _frame;
    char _gene_id[GTF_MAX_READ_SIZE];
    char _transcript_id[GTF_MAX_READ_SIZE];
    
    //Tried stringstream - it is horrible, not thread safe, incredibly slow in a multithreaded application
    sscanf(line.c_str(), "%s\t%s\t%s\t%u\t%u\t%s\t%c\t%c\t%*[^\"]\"%255[^\"]\"%*[^\"]\"%255[^\"]\"", _chr, _source, _feature, &_start, &_end, _score, &_strand, &_frame, &_gene_id, &_transcript_id);
        
    chr = _chr;
    source = _source;
    feature = _feature;
    start = _start;
    end = _end;
    score = _score;
    strand = _strand;
    frame = _frame;
    gene_id = _gene_id;
    transcript_id = _transcript_id;
    
}
*/
