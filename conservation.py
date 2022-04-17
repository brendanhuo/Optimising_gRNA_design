from utils import *
from Bio.Blast import NCBIWWW, NCBIXML

def runBlast(seq, database = "nt", numbers = 6000):
    """returns numbers amount of results. Issues where if ask for too many results then XML output breaks and can't be parsed so reduce number if see issue"""
    #databases available - https://ncbi.github.io/blast-cloud/blastdb/available-blastdbs.html
    result_handle = NCBIWWW.qblast("blastn", database, seq, expect = 1000, descriptions = numbers/2, alignments = numbers/2, hitlist_size = numbers, short_query = True)
    return result_handle

def saveBlasttoXML(result_handle, file_location):
    with open(file_location, 'w') as save_file: 
      blast_results = result_handle.read() 
      save_file.write(blast_results)

def runBlastforallgRNAS(records,file_location, search_for = 'Hepatitis B',minNumber = 1500):
    for i in range(len(records)):
        save_file_location = file_location + "/gRNA%s.xml" %i
        total_length = 0
        print('gRNA%s' %i)
        try:
            for record in NCBIXML.parse(open(save_file_location)): 
                for align in record.alignments:
                    if search_for in align.title and 'isolate' in align.title:
                        for hsp in align.hsps:
                            total_length += hsp.align_length
            if total_length == 0:
                print("gRNA%s - fetching again, 0 results issue" % i)
                result_handle = runBlast(records[i].seq, numbers = minNumber//2)
                saveBlasttoXML(result_handle, save_file_location)
                print("gRNA%s done" % i)
            else:
                print('gRNA check okay')

        except:
            print("gRNA%s - fetching" % i)
            result_handle = runBlast(records[i].seq, numbers = minNumber)
            saveBlasttoXML(result_handle, save_file_location)
            print("gRNA%s done" % i)

            for record in NCBIXML.parse(open(save_file_location)): 
                for align in record.alignments:
                    if search_for in align.title and 'isolate' in align.title:
                        for hsp in align.hsps:
                            total_length += hsp.align_length
            if total_length == 0:
                print("gRNA%s - fetching again, 0 results issue" % i)
                result_handle = runBlast(records[i].seq, numbers = minNumber//2)
                saveBlasttoXML(result_handle, save_file_location)
                print("gRNA%s done" % i)

def computeConservation(records,file_location,search_for = 'Hepatitis B', plot = True, save = True):
    conservation_scores = []
    for i in range(len(records)):
        save_file_location = file_location + "/gRNA%s.xml" %i
        count = 0
        total_correct = 0
        total_length = 0
        try:
            for record in NCBIXML.parse(open(save_file_location)): 
                for align in record.alignments:
                    if search_for in align.title and 'isolate' in align.title:
                        for hsp in align.hsps:
                            count += 1
                            total_correct += hsp.positives
                            total_length += hsp.align_length
                    else:
                        pass
            print("gRNA%s" % i)
            conservation_scores.append(total_correct/total_length)
        except:
            print("gRNA%s - rip didn\'t work" % i)
            conservation_scores.append(None)
    
    if save:
        conservation_scores_np = np.asarray(conservation_scores)
        cons_save_file_location = file_location + "/conservation_scores.npy"
        np.save(cons_save_file_location,conservation_scores_np)
    
    if plot:
        for i in range(len(conservation_scores)):
            if conservation_scores[i] == None:
                conservation_scores[i] = 0
        x = np.arange(len(conservation_scores))
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        ax.bar(x, conservation_scores, width = 0.6, color = 'r')
        plt.xlabel('gRNA sequence');plt.ylabel('Average % conserved in Hep B virus genome');plt.grid()
        plt.ylim(0.92, 1), plt.xlim(0, len(records))
        plt.show()
    
    return conservation_scores_np