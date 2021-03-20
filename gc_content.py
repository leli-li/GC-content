# -*- coding: utf-8 -*-
"""
This program is a tool to plot GC content
@author: lily
"""

#import standard modules
import sys
import os
import argparse

#import non-standard modules
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

#check if the non-standard modules are installed
try:
    import matplotlib
except ModuleNotFoundError:
    print('No module named matplotlib'+'\n'+'please install it')
    sys.exit()


# create the parser

gc_parser = argparse.ArgumentParser(description='A tool for plotting GC content based on sliding window')

#add position argument: mode argument
#this program has five modes of processing files
gc_parser.add_argument('mode',help='choose mode of program: single(s), pairwise comparison(c), batch processing(b), batch pairwise comparison(bc) and all lines in one plot(a)', choices=["s","c","b","bc","a"],nargs=1,type=str)

#add other arguments
gc_parser.add_argument("-f","--file", help="input file", nargs="+", type=str) #
gc_parser.add_argument("-w","--window", help="set the size of sliding window", nargs=1, type=int)
gc_parser.add_argument("-s","--step", help="set the step", nargs=1, type=int)
gc_parser.add_argument("-r","--range",metavar=("Lower limit","upper limit"),help="set the range of GC content to filter sequences. between 0 and 100", nargs=2, type=int)
# Execute the parse_args() method
args = gc_parser.parse_args()


def file_processing():
#this function is to store the id and corresponding sequence into a dictionary
#this dictionary is used to plot the distribution of gc content of each sequence 
#and used to filter sequences by a range of gc content 
#also combine sequences to a whole one. this is used to plot gc content through the whole sequence(usually the genome)
    global seq_dict     #make variables to be used in the whole program
    global f_basename
    global whole_seq
    f = open(file,"r") 
    f_basename=file.split(".")[0] #extract file name without suffix
    seq_dict={} #the dictionary of id and corresponding sequence
    for line in f:  #combine the sequence if the sequence occupies more than one line
        if ">" in line:
            seq_id = line.split(" ")[0]
            seq_id = seq_id.strip()
            seq_id = seq_id.strip(">")
            seq = ""               
        if ">" not in line:
            line=line.strip()
            line=line.strip("-")
            line=line.strip("-")
            seq=seq+line           
            seq_dict[seq_id]=seq  #key is the id and value is the sequence
#print some information to tell the user about running process
    print("\n" + "File:" + f_basename+"\n")  
    print("Number of sequences:"+ str(len(seq_dict))+"\n")
#combine sequence to a whole sequence  
#it is used to plot GC content through the whole sequence(usually the genome)  
    whole_seq=""
    for seq in seq_dict.values():
        whole_seq=whole_seq+seq
    whole_seq=whole_seq.upper()
    print("Total length:"+ str(len(whole_seq))) 

def GC_content_distribution():
#this function is to calculate GC content of each sequence and plot histogram:
    seq_gc=[]  #store GC contents of each sequence 
    for id in seq_dict:
        numerator=seq_dict[id].upper().count('G')+seq_dict[id].upper().count('C') #when calculating GC content, convert letters to uppercase letters
        denominator=seq_dict[id].upper().count('G')+seq_dict[id].upper().count('C')+seq_dict[id].upper().count('A')+seq_dict[id].upper().count('T')
        if denominator!=0:
            ratio = numerator/denominator
            seq_gc.append(ratio)
        if denominator==0: #consider the situation of there's no stand nucleotides in a sequence
            ratio = 0
            seq_gc.append(0)
    fig = plt.figure(figsize=(15,10),dpi=300) #set the size of graph , it is only a blank graph
    ax = fig.add_subplot(1, 1, 1) #only add one figure
    ax.hist(seq_gc, bins=10, histtype='bar',color='steelblue') #plot
    plt.xlim(0,1) #set the range of x axis
    plt.xlabel('GC content',fontsize = 20) #add the labels of axis
    plt.ylabel('Counts of sequences',fontsize = 20)
    plt.tick_params(labelsize=15) #set the size of valus in axis
    plt.title(f_basename , fontsize=30,fontweight="heavy") #set the title of graph
    plt.savefig(f_basename+"_fre.png")   #save graph , the name is based on input file name



def calculate_GC():
#this function is to calculate GC content based on sliding window
    global window_size   #make variables to be used in the whole program
    global step
    global window
    global remainder
    global middle_p
    global interval
    global gc_list
    global stop_list
    global remainder_gc
    window_size = int(args.window[0]) #retrieve the window size from command
    step = int(args.step[0])
    #set the value of interval of x axis
    interval = str(int(len(whole_seq)/6))
    inter_list=list(interval)
    for i in range(1,len(inter_list)): #make the other digits are zero except the first one, so that the value in x axis are the numbers ending with zero
        inter_list[i]="0"
    interval=int("".join(inter_list))    
    seq_length = len(whole_seq)
    window={}
    remainder=[]    
    stop_list=[]
    middle_p=[] 
    start_list=[]
    if step != 0: #consider different situations of step value
        for start_index in range(0,seq_length,step): #the step decides the start position of a sliding window ,and the window size decides the length ,that is the stop position of a sliding window
            stop_index = start_index + window_size -1
            if stop_index <= seq_length-1: 
                window[start_index] = stop_index
                stop_list.append(stop_index)  
                start_list.append(start_index)
    if step ==0:
        for start_index in range(0,seq_length,window_size):            
            stop_index = start_index + window_size -1
            if stop_index <= seq_length-1: 
                window[start_index] = stop_index
                stop_list.append(stop_index)   
    #deal with the remainder part in a sequence: 
    if step != 0:
        if stop_list[-1] < seq_length -1:  
            for a in range(start_list[-1]+step, seq_length):
                remainder.append(a)
    if step == 0:
        if stop_list[-1] < seq_length -1 :
            for a in range(int(stop_list[-1])+1, seq_length):
                remainder.append(a)
#generate the values in x axis, which is the middle position of a window
#the index should plus 1 because python index starts with 0, but in real sequence index starts with 1
    for index in window.keys():
        middle = ((index+1)+(window[index]+1))/2 
        middle_p.append(middle) 
    if len(remainder)>0:
        remainder_middle=((remainder[0]+1)+(remainder[-1]+1))/2
        middle_p.append(remainder_middle)        
#calculate_gc_content:    
    gc_list = []  #this list is the y axis 
    for a,b in zip (list(window.keys()),list(window.values())):
        segment=whole_seq[a:(b+1)]
        if (segment.count("G")+segment.count("C")+segment.count("A")+segment.count("T")) !=0:
            gc=(segment.count("G")+segment.count("C")) / (segment.count("G")+segment.count("C")+segment.count("A")+segment.count("T"))
            gc_list.append(gc)
        if (segment.count("G")+segment.count("C")+segment.count("A")+segment.count("T")) ==0:
            gc_list.append(0)
    if len(remainder)>0:
        remainder_segment = whole_seq[remainder[0]:]
        if (remainder_segment.count("G")+remainder_segment.count("C")+remainder_segment.count("A")+remainder_segment.count("T")) !=0:
            remainder_gc = (remainder_segment.count("G")+remainder_segment.count("C")) / (remainder_segment.count("G")+remainder_segment.count("C")+remainder_segment.count("A")+remainder_segment.count("T")) 
            gc_list.append(remainder_gc)
        if (remainder_segment.count("G")+remainder_segment.count("C")+remainder_segment.count("A")+remainder_segment.count("T")) ==0:
            gc_list.append(0)

def single_line_chart():
#plot of one single file 
#generate a baseline of mean value
#calculate mean value of sliding windows
    total_gc=0
    for value in gc_list:
        total_gc=total_gc+value
    mean_gc=total_gc/len(gc_list)        
    x = middle_p
    y = gc_list
    fig = plt.figure(figsize=(20,10),dpi=300) #set a blank figure
    ax = fig.add_subplot(1, 1, 1) # add one plot
    ax.plot(x,y,label=f_basename,linewidth=2.3) #plot the line of GC content
    ax.xaxis.set_major_locator(ticker.MultipleLocator(interval)) #set the values displayed of x axis
    ax.ticklabel_format(useOffset=False, style='plain') #not show the scientific notation
    #plot a straight line of mean value
    a = [0,len(whole_seq)] 
    b = [mean_gc,mean_gc]
    ax.plot(a,b,label="mean="+str(round(mean_gc,2)),color="lightcoral",linewidth=2,linestyle="--")
    plt.title(f_basename,fontsize=30,fontweight="heavy") #set the title
    plt.xlabel('Position',fontsize = 20) #set the label of axis
    plt.ylabel('GC content',fontsize = 20)
    plt.ylim(min(gc_list)-0.1,max(gc_list)+0.1)
    plt.tick_params(labelsize=15) #set the size of numbers of axis
    plt.legend(ncol=1,fontsize=15,loc="upper right") #set the size of legend and make it show in one column
    plt.text(len(whole_seq)/2,min(gc_list)-0.1+0.02,s="window size="+str(window_size)+"\n"+"step="+str(step),horizontalalignment='center',verticalalignment='top',fontsize=12,color="slategrey") #add the text of window size and step
    plt.savefig(f_basename+"_line.png") #save the graph



def sequence_filter():
#this function is to filter sequence of a specific range  
#get the range from the arguments in command
    lower_limit=float(args.range[0])
    upper_limit=float(args.range[1])
#calculate GC content of each sequence and write to file
    out=open(f_basename+"_filter.out","w")
    for id in seq_dict:
        numerator=seq_dict[id].upper().count('G')+seq_dict[id].upper().count('C')
        denominator=seq_dict[id].upper().count('G')+seq_dict[id].upper().count('C')+seq_dict[id].upper().count('A')+seq_dict[id].upper().count('T')
        if denominator!=0:
            ratio = round(numerator/denominator,2)
            if lower_limit <= ratio *100 <= upper_limit:
                out.write(">"+id+"\t"+"GC="+str((ratio))+"\n"+seq_dict[id]+"\n")
        if denominator==0:
            ratio = 0
            if lower_limit <= ratio * 100 <= upper_limit:
                out.write(">"+id+"\t"+"GC="+str((ratio))+"\n"+seq_dict[id]+"\n")
    out.close()   


def another_file_processing():
#this function is used in pairwise comparison
#it is used to process another file
    global seq_dict2
    global f_basename2
    global whole_seq2
    f = open(file2,"r")    
    f_basename2=file2.split(".")[0]
    seq_dict2={}
    for line in f:  #combine the sequence if the sequence occupies more than one line
        if ">" in line:
            seq_id = line.split(" ")[0]
            seq_id = seq_id.strip()
            seq_id = seq_id.strip(">")
            seq = ""   
        if ">" not in line:
            line=line.strip()
            line=line.upper()
            seq=seq+line
            seq_dict2[seq_id]=seq 
#print some information to tell the user about running process
    print("\n" + "File2:" + f_basename2+"\n")
    print("Number of sequences:"+ str(len(seq_dict2))+"\n")
    f.close()
#combine sequence to a whole sequence    
    whole_seq2=""
    for seq in seq_dict2.values():
        whole_seq2=whole_seq2+seq
    whole_seq2=whole_seq2.upper()
    print("Total length:"+ str(len(whole_seq2)))


def another_GC_and_comparison_plot():
#this function is used in pairwise comparison
#calculate GC_content of file2 and plot line charts of file1 and file2 in a same figure
    global window2
    global remainder2
    global middle_p2
    global interval2
    global gc_list2
    #calculate the interval of x axis
    interval2 = str(int(len(whole_seq2)/6))
    inter_list=list(interval2)
    for i in range(1,len(inter_list)): #make the other digits are zero except the first one, so that the value in x axis are the numbers ending with zero
        inter_list[i]="0"
    interval2=int("".join(inter_list))    
    seq_length = len(whole_seq2)
    window2={}
    remainder2=[]    
    stop_list=[]
    middle_p2=[] 
    start_list=[]
    if step != 0:
        for start_index in range(0,seq_length,step):
            stop_index = start_index + window_size -1
            if stop_index <= seq_length-1: 
                window2[start_index] = stop_index
                stop_list.append(stop_index)
                start_list.append(start_index)
    if step ==0:
        for start_index in range(0,seq_length,window_size):            
            stop_index = start_index + window_size -1
            if stop_index <= seq_length-1: 
                window2[start_index] = stop_index
                stop_list.append(stop_index)    
#deal with the remainder part in sequence:    
    if step != 0:
        if stop_list[-1] < seq_length -1 : 
            for a in range(start_list[-1]+step, seq_length):
                remainder2.append(a)
    if step == 0:
        if stop_list[-1] < seq_length -1 :
            for a in range(int(stop_list[-1])+1, seq_length):
                remainder2.append(a)
#generate the x axis, which is the middle position of a window
#the index should plus 1 because python index starts with 0, but in real sequence index starts with 1
    for index in window2.keys():
        middle = ((index+1)+(window2[index]+1))/2 
        middle_p2.append(middle) 
    if len(remainder2)>0:
        remainder_middle=((remainder2[0]+1)+(remainder2[-1]+1))/2
        middle_p2.append(remainder_middle)        
#calculate_gc_content:    
    gc_list2 = []
    for a,b in zip (list(window2.keys()),list(window2.values())):
        segment=whole_seq2[a:(b+1)]
        if (segment.count("G")+segment.count("C")+segment.count("A")+segment.count("T")) !=0:
            gc=(segment.count("G")+segment.count("C")) / (segment.count("G")+segment.count("C")+segment.count("A")+segment.count("T"))
            gc_list2.append(gc)
        if (segment.count("G")+segment.count("C")+segment.count("A")+segment.count("T")) ==0:
            gc_list2.append(0)
    if len(remainder2)>0:
        remainder_segment = whole_seq2[remainder2[0]:]
        if (remainder_segment.count("G")+remainder_segment.count("C")+remainder_segment.count("A")+remainder_segment.count("T")) !=0:
            remainder_gc = (remainder_segment.count("G")+remainder_segment.count("C")) / (remainder_segment.count("G")+remainder_segment.count("C")+remainder_segment.count("A")+remainder_segment.count("T")) 
            gc_list2.append(remainder_gc)
        if (remainder_segment.count("G")+remainder_segment.count("C")+remainder_segment.count("A")+remainder_segment.count("T")) ==0:
            gc_list2.append(0)

#calculate average GC content of genome1
    total_gc=0
    for value in gc_list:
        total_gc=total_gc+value
    mean_gc=total_gc/len(gc_list)
#calculate average GC content of genome2
    total_gc2=0
    for value in gc_list2:
        total_gc2=total_gc2+value
    mean_gc2=total_gc2/len(gc_list2)
#generate a baseline: average GC content of the two mean values
    ave_gc=round((mean_gc+mean_gc2)/2,2)  
    fig = plt.figure(figsize=(20,10),dpi=300)
    ax = fig.add_subplot(1, 1, 1)
#plot the first genome        
    x = middle_p
    y = gc_list   
    ax.plot(x,y,label=f_basename+":mean="+str(round(mean_gc,2)),linewidth=2.3)   
    ax.ticklabel_format(useOffset=False, style='plain') #not show the scientific notation
#plot the second genome 
    x2 = middle_p2
    y2 = gc_list2 
    ax.plot(x2,y2,label=f_basename2+":mean="+str(round(mean_gc2,2)),linewidth=2.3)
    ax.ticklabel_format(useOffset=False, style='plain') #not show the scientific notation
#plot the baseline
    a = [0,max(len(whole_seq),len(whole_seq2))]
    b = [ave_gc,ave_gc]
    ax.plot(a,b,label="mean="+str(round(ave_gc,2)),color="lightcoral",linewidth=2,linestyle="--")
#set the parameters of the plot
    ax.xaxis.set_major_locator(ticker.MultipleLocator(max(interval,interval2))) #the interval of x axis should be the max value
    plt.ylim(min(min(gc_list),min(gc_list2))-0.1,max(max(gc_list),max(gc_list2))+0.1)
    plt.title(f_basename+" VS "+f_basename2,fontsize=30,fontweight="heavy")
    plt.xlabel('Position',fontsize = 20)
    plt.ylabel('GC content',fontsize = 20)
    plt.tick_params(labelsize=15)
    plt.legend(ncol=1,fontsize=15,loc="upper right")
    #add the text of window size and step, and make its position in the middle bottom of the graph
    plt.text(max(len(whole_seq),len(whole_seq2))/2,min(min(gc_list),min(gc_list2))-0.1+0.02,s="window size="+str(window_size)+"\n"+"step="+str(step),horizontalalignment='center',verticalalignment='top',fontsize=12,color="slategrey")
    plt.savefig(f_basename+"_VS_"+f_basename2+".png")

#The following part is to use the function with different arguments
#check if the file is missing
if args.mode[0]== "s" and not args.file:   
    print("input file missing, nothing will happen")
    sys.exit()


#mode of single file
#input file and processing
if args.mode[0]=="s" and args.file: #if this argument is used
    file = args.file[0]
    try:                        #check the existence of file
        open(file,"r")
    except IOError:
        print("file does not exit,please check the file name")
        sys.exit()   
    file_processing() 
    print("\n"+"File is read successfully!"+"\n")
    if args.range: #use the range argument to filter sequences
        print("sequence filtering is processing" + "\n")
        if args.range[0]> args.range[1]:
            print ("lower limit cannot be higher than upper limit!")
            sys.exit()
        else:
            sequence_filter()
            print("Output file is saved as " + f_basename+"_filter.out")
    if not args.range: 
        if not args.window or not args.step:
            print( "Histogram of GC content distribution will be plotted automatically" + "\n")
            GC_content_distribution()
            print("Histogram is saved as " + f_basename +"_dist.png")
        if args.window and args.step:  #use the argument of window size and step
            print("Line chart is being processed" +"\n")
            calculate_GC()
            single_line_chart()
            print("Line chart is saved as " + f_basename +"_line.png" + "\n")
    
            
                      
#mode of pairwise comparison  
#check if the file is missing
if args.mode[0]=="c" and not args.file:   
    print("input file missing, nothing will happen")
    sys.exit()
    
if args.mode[0]=="c" and len(args.file) ==1:   
    print("Two files are needed, please input another file")
    sys.exit()
    
if args.mode[0]=="c" and args.file: #if this argument is used
    file = args.file[0]
    try:                        #check the existence of file
        open(file,"r")
    except IOError:
        print("file does not exit,please check the file name")
        sys.exit()   
    file_processing()
    file2 = args.file[1]
    try:                        #check the existence of file
        open(file2,"r")
    except IOError:
        print("file2 does not exit,please check the file name")
        sys.exit()   
    another_file_processing()
    if args.window and args.step:
        print("Line chart is being processed")
        calculate_GC()
        another_GC_and_comparison_plot()
        print("Line chart is saved as " + f_basename+"_VS_"+f_basename2 +".png" + "\n")
        
        
#mode of batch processing
#this mode is to generate line charts of multiple files
if args.mode[0]=="b":
#retrive the path of this script
#when running the script, it should be put in the same directory of the target files
    script_path=os.path.dirname(os.path.realpath('__file__'))
    target_file=[] #index of files to be processed
    for root,dirs,files in os.walk(script_path):
        for file in files:
            if "py" not in file:
                target_file.append(file)
    for file in target_file:
        file_processing()
        print("\n"+"File is read successfully!"+"\n")
        if args.window and args.step:
            print("Line chart is being processed" +"\n")
            calculate_GC()
            single_line_chart()
            print("Line chart is saved as " + f_basename +"_line.png" + "\n")
        

#mode of batch comparison
#this mode is batch processing of pairwise comparison          
if args.mode[0]=="bc":
#retrive the path of the script
#when running the script, it should be put in the same directory of the target files
    script_path=os.path.dirname(os.path.realpath('__file__'))
    target_file=[] #index of files to be processed
    for root,dirs,files in os.walk(script_path):
        for file in files:
            if "py" not in file:
                target_file.append(file)    
    if not args.window or not args.step:
        print("please set the window size and step")
    if args.window and args.step:
        for i in range(0,len(target_file)):
            for j in range(i+1,len(target_file)):
                file = target_file[i]
                file2 = target_file[j]
                file_processing()
                another_file_processing()
                calculate_GC()
                another_GC_and_comparison_plot()
                print("Line chart is saved as " + f_basename+"_VS_"+f_basename2 +".png" + "\n")
                
            

#the following part is to plot all lines in one figure   
def all_in_one_plot():
#this function is to draw a simple line of each file
    global mean_gc
    total_gc=0
    for value in gc_list:
        total_gc=total_gc+value
    mean_gc=total_gc/len(gc_list)        
    x = middle_p
    y = gc_list
    ax.plot(x,y,label=f_basename+":mean="+str(round(mean_gc,2)),linewidth=2.3)
            
if args.mode[0]=="a": 
    script_path=os.path.dirname(os.path.realpath('__file__'))
    target_file=[] #index of files to be processed
    for root,dirs,files in os.walk(script_path):
        for file in files:
            if "py" not in file:
                target_file.append(file)          
    if not args.window or not args.step:
        print("please set the window size and step")
    if args.window and args.step:
        fig = plt.figure(figsize=(20,10),dpi=300) #generate a blank  figure , add a line every time the file is read
        ax = fig.add_subplot(1, 1, 1)
        total_mean = 0
        interval_assembly=[]
        whole_length_assembly=[]
        min_gc=[]
        max_gc=[]
        for file in target_file:
            file_processing() 
            calculate_GC()
            all_in_one_plot()
            total_mean = total_mean + mean_gc 
            interval_assembly.append(interval)  #store important values in list
            whole_length_assembly.append(len(whole_seq))
            min_gc.append(min(gc_list)) 
            max_gc.append(max(gc_list))
        total_mean = total_mean / len(target_file)
        #plot baseline of mean value
        a = [0,max(whole_length_assembly)]
        b = [total_mean,total_mean]
        ax.plot(a,b,label="total_mean="+str(round(total_mean,2)),color="lightcoral",linewidth=2,linestyle="--")
        ax.xaxis.set_major_locator(ticker.MultipleLocator(max(interval_assembly))) #the interval value is the max one
        ax.ticklabel_format(useOffset=False, style='plain') #not show the scientific notation
        plt.ylim(min(min_gc)-0.1,max(max_gc)+0.1)
        plt.xlabel('Position',fontsize = 20)
        plt.ylabel('GC content',fontsize = 20)
        plt.tick_params(labelsize=15)
        plt.legend(ncol=1,fontsize=12,loc="upper right") 
        #add the text and make its position in the middle bottom of the graph
        plt.text(max(whole_length_assembly)/2,min(min_gc)-0.1+0.02,s="window size="+str(window_size)+"\n"+"step="+str(step),horizontalalignment='center',verticalalignment='top',fontsize=12,color="slategrey")
        plt.savefig("all_line.png") 
        print("Line chart is saved as all_line.png")
    

