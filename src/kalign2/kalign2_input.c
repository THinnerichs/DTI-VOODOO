/*
	kalign2_input.c 
	
	Released under GPL - see the 'COPYING' file   
	
	Copyright (C) 2006 Timo Lassmann <timolassmann@gmail.com>
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    
	Please send bug reports, comments etc. to:
	timolassmann@gmail.com
*/

#include "kalign2.h"
#include "kalign2_input.h"

#include <ctype.h>


struct alignment* detect_and_read_sequences(struct alignment* aln,struct parameters* param)
{
	
	int feature = 0;
	char **input = 0;
	unsigned short int* input_type = 0;
	unsigned short int* input_numseq = 0;
	
	int num_input = 0;
	int i = 0;
	int j = 0;
	int c = 0;
	int free_read = 1;
	while(free_read == 1 || param->infile[i]){
		num_input++;
		i++;
		free_read = 0;
	}
	numseq = 0;

	
	input = malloc(sizeof(char*) * num_input);
	input_type = malloc(sizeof(unsigned short int) * num_input);
	input_numseq = malloc(sizeof(unsigned short int) * num_input);
	
	for (i = 0; i < num_input;i++){
		input[i] = 0;
		input_type[i] = 0;
		input_numseq[i] = 0;
	}

	free_read = 0;
	
	if(param->quiet){
		c = 1;
	}else{
		c = 0;
	}
	
	
	for (i = c; i < num_input;i++){
		if(!param->infile[i]){
			fprintf(stderr,"reading from STDIN: ");
		}else{
			fprintf(stderr,"reading from %s: ",param->infile[i]);
		}
		input[i] = get_input_into_string(input[i],param->infile[i]);
		if(input[i]){
			free_read++;
			if (byg_start("<macsim>",input[i]) != -1){
				input_numseq[i] = count_sequences_macsim(input[i]);
				feature = 1;
				input_type[i] = 1;
			}else if (byg_start("<uniprot",input[i]) != -1){
				input_numseq[i] = count_sequences_uniprot(input[i]);
				input_type[i] = 2;
			}else if(byg_start("This SWISS-PROT",input[i]) != -1){
				input_numseq[i] = count_sequences_swissprot(input[i]);
				input_type[i] = 3;
			}else if (byg_start("This Swiss-Prot",input[i]) != -1){
				input_numseq[i] = count_sequences_swissprot(input[i]);
				input_type[i] = 3;
			}else if (byg_start("CLUSTAL W",input[i]) != -1){
				input_numseq[i] = count_sequences_clustalw(input[i]);
				input_type[i] = 4;
			}else if (byg_start("PileUp",input[i]) != -1){
				input_numseq[i] = count_sequences_clustalw(input[i]);
				input_type[i] = 4;
			}else if (byg_start("MSF:",input[i]) != -1){
				input_numseq[i] = count_sequences_clustalw(input[i]);
				input_type[i] = 4;
			}else if (byg_start("STOCKHOLM",input[i]) != -1){
				input_numseq[i] = count_sequences_stockholm(input[i]);
				input_type[i] = 5;
			}else{
				input_numseq[i]  = count_sequences_fasta(input[i]);
				input_type[i] = 0;
			}
			fprintf(stderr,"found %d sequences\n",input_numseq[i]);
			
			if(input_numseq[i] < 1){
				free(input[i]);
				input[i] = 0;
			}else{
				numseq += input_numseq[i];
			}
		}else{
			fprintf(stderr,"found no sequences.\n");
			if(!param->outfile && i){
				param->outfile = param->infile[i];
				fprintf(stderr,"-> output file, in ");
				//try to set format.... 
				if(!param->format){
					if (byg_start("msf",param->outfile) != -1){
						param->format = "msf";
					}else if (byg_start("clustal",param->outfile) != -1){
						param->format = "clustal";
					}else if (byg_start("aln",param->outfile) != -1){
						param->format = "clustal";
					}else if (byg_start("macsim",param->outfile) != -1){
						param->format = "macsim";
					}else{
						param->format = "fasta";
					}
					if(param->reformat){
						fprintf(stderr,"unaligned fasta format\n");
					}else if(param->format){
						fprintf(stderr,"%s format\n",param->format);
					}else{
						fprintf(stderr,"fasta format\n");
					}
				}
			}
			fprintf(stderr,"\n");
		}
	}

	
	if(numseq < 2){
		fprintf(stderr,"%s\n", usage);
		if(!numseq){
		fprintf(stderr,"\nWARNING: No sequences found.\n\n");
		}else{
		fprintf(stderr,"\nWARNING: Only one sequence found.\n\n");
		}
		for (i = 0; i < num_input;i++){
			free(input[i]);
		}
		free(input_numseq);
		free(input_type);
		free(input);
		free_param(param);
		exit(0);
	}

	if(byg_start(param->alignment_type,"profPROFprofilePROFILE") != -1){
		if( free_read  < 2){
			fprintf(stderr,"\nWARNING: You are trying to perform a profile - profile alignment but ony one input file was detected.\n\n");
			param->alignment_type = "default";
		}
	}

	
	if (param->feature_type && !feature){
		fprintf(stderr,"\nWARNING: You are trying to perform a feature alignment but the input format(s) do not contain feature information.\n");
		for (i = 0; i < num_input;i++){
			free(input[i]);
		}
		free(input_numseq);
		free(input_type);
		free(input);
		free_param(param);
		exit(0);
	}
	
	numprofiles = (numseq << 1) - 1;
	aln = aln_alloc(aln);
	//numseq = 0;
	if(byg_start(param->alignment_type,"profPROFprofilePROFILE") != -1){
		j = 0;
		for (i = 0; i < num_input;i++){
			
			if(input[i]){
					
				switch(input_type[i]){
					case 0:
						aln = read_alignment(aln,input[i]);
						break;
					case 1:
						aln = read_alignment_macsim_xml(aln,input[i]);
						break;
					case 2:
						aln = read_alignment_uniprot_xml(aln,input[i]);
						break;
					case 3:

						aln = read_alignment_from_swissprot(aln, input[i]);
						break;
					case 4:
						aln = read_alignment_clustal(aln,input[i]);
						break;
					case 5:
						aln = read_alignment_stockholm(aln,input[i]);
						break;
					
					default:
						aln = read_alignment(aln,input[i]);
						break;
				}
				input[i] = 0;
				//create partial profile....
				aln->nsip[numseq+j] = input_numseq[i];
				aln->sip[numseq+j] = malloc(sizeof(int)*aln->nsip[numseq+j]);
				
				//fprintf(stderr,"%d	%d\n",numseq+j,aln->sl[numseq+j]);
				j++;
			}
		}
		num_input = j;
		c = 0;
		for (i = 0;i < num_input;i++){
		//	
			for ( j = 0; j < aln->nsip[numseq+i];j++){
				aln->sip[numseq+i][j] = c;
				c++;
		//		fprintf(stderr,"%d ",aln->sip[numseq+i][j]);
			}
			aln->sl[numseq+i] = aln->sl[aln->sip[numseq+i][0]];
		//	fprintf(stderr,"PROFILE:%d	contains: %d long:%d\n",i+numseq,aln->nsip[numseq+i],aln->sl[numseq+i]);
	//		fprintf(stderr,"\n");
		}
		
		//sanity check -are all input 
		int a,b;
		for (i = 0;i < num_input;i++){
			for ( j = 0; j < aln->nsip[numseq+i]-1;j++){
				a = aln->sip[numseq+i][j];
				a = aln->sl[a];
				for (c =  j+1; j < aln->nsip[numseq+i];j++){
					b = aln->sip[numseq+i][c];
					b = aln->sl[b];
					if(a != b){
						fprintf(stderr,"Unaligned sequences in input %s.\n",param->infile[i]);
						for (i = 0; i < num_input;i++){
							free(input[i]);
						}
						free(input_numseq);
						free(input_type);
						free(input);
						free_aln(aln);
						free_param(param);
						exit(0);
					}
				}
				
			}

		}
		
		//exit(0);
		
		/*for (i = 0; i < numseq;i++){
			fprintf(stderr,"len%d:%d\n",i,aln->sl[i]);	
			for ( j =0 ; j < aln->sl[i];j++){
				//if(aln->s[i][j]> 23 || aln->s[i][j] < 0){
				//	 aln->s[i][j] = -1;
				//}
				fprintf(stderr,"%d ",aln->s[i][j]);
			}
		//	fprintf(stderr,"\n");
		}
		exit(0);*/
	}else{
		for (i = 0; i < num_input;i++){
			if(input[i]){
				switch(input_type[i]){
					case 0:
						aln = read_sequences(aln,input[i]);
						break;
					case 1:
						aln = read_sequences_macsim_xml(aln,input[i]);
						break;
					case 2:
						aln = read_sequences_uniprot_xml(aln,input[i]);
						break;
					case 3:
						aln = read_sequences_from_swissprot(aln, input[i]);
						break;
					case 4:
						aln = read_sequences_clustal(aln,input[i]);
						break;
					case 5:
						aln = read_sequences_stockholm(aln,input[i]);
						break;
					
					default:
						aln = read_sequences(aln,input[i]);
						break;
				}
				/*if (byg_start("<macsim>",input[i]) != -1){
					aln = read_sequences_macsim_xml(aln,input[i]);
				}else if (byg_start("<uniprot",input[i]) != -1){
					aln = read_sequences_uniprot_xml(aln,input[i]);
				}else if(byg_start("This SWISS-PROT entry is copyright.",input[i]) != -1){
					aln = read_sequences_from_swissprot(aln, input[i]);
				}else if (byg_start("This Swiss-Prot entry is copyright.",input[i]) != -1){
					aln = read_sequences_from_swissprot(aln, input[i]);
				}else if (byg_start("CLUSTAL W",input[i]) != -1){
					aln = read_sequences_clustal(aln,input[i]);
				}else if (byg_start("PileUp",input[i]) != -1){
					aln = read_sequences_clustal(aln,input[i]);
				}else if (byg_start("MSF:",input[i]) != -1){
					aln = read_sequences_clustal(aln,input[i]);
				}else if (byg_start("STOCKHOLM",input[i]) != -1){
					aln = read_sequences_stockholm(aln,input[i]);
				}else{
					aln = read_sequences(aln,input[i]);
				}*/
				input[i] = 0;
			}
		}
	}
	if(numseq < 2){
		fprintf(stderr,"\nNo sequences could be read.\n");
		free_param(param);
		exit(0);
	}
	if(!param->format && param->outfile){
			if (byg_start("msf",param->outfile) != -1){
				param->format = "msf";
			}else if (byg_start("clustal",param->outfile) != -1){
				param->format = "clustal";
			}else if (byg_start("aln",param->outfile) != -1){
				param->format = "clustal";
			}else if (byg_start("macsim",param->outfile) != -1){
				param->format = "macsim";
			}
			fprintf(stderr,"Output file: %s, in %s format.\n",param->outfile,param->format);
	}
	
	
	free(input);
	free(input_type);
	free(input_numseq);
	return aln;
}

int count_sequences_macsim(char* string)
{
	int n = 0;
	n = byg_count("<seq-name>",string);
	if(!n){
		return -1;
	}
	return n;
}

int count_sequences_swissprot(char* string)
{
	int n = 0;
	n = byg_count("ID   ",string);
	if(!n){
		return 0;
	}
	return n;
}

int count_sequences_uniprot(char* string)
{
	int n = 0;
	n = byg_count("<entry",string);
	if(!n){
		return 0;
	}
	return n;
}

int count_sequences_stockholm(char* string)
{
	char* p1 = string;
	int i = 0;
	int j = 0;
	int n = 0;
	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		if (!(byg_start("//",p1))){
			break;
		}	
		j = byg_end("#",p1);
		if(j != 1){
			n++;
		}
	}
	if(!n){
		return 0;
	}
	return n;
}

int count_sequences_clustalw(char* string)
{
	char* p1 = string;
	int i = 0;
	int j = 0;
	int c = 0;
	int n = 0;
	int f = 0;
	

	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		j = byg_end(" ",p1);
		f = byg_end("\n",p1);
		if(f > 2 && f>j && j!= 1){
			if(c ==0){
				i = j;
				while(p1[i] != '\n'){
					//if (!isspace((int)p1[i])){
					//	len++;
					//}
					i++;
				}		
			}
			c++;
		}else{
			if (c){
				if(c > n){
					n = c;
				}
				c =0;
			}
		}
	}
	if(!n){
		return 0;
	}
	return n;
}

int count_sequences_fasta(char* string)
{
	int nbytes;
	int i; 
	int n = 0;
	int stop = 0;
	nbytes = strlen(string);
	for (i =0;i < nbytes;i++){
		if (string[i] == '>'&& stop == 0){
			stop = 1;
			n++;
		}
		if (string[i] == '\n'){
			stop = 0;
		}
	}
	if(!n){
		return 0;
	}
	return n;
}




char* get_input_into_string(char* string,char* infile)
{
	int i = 0;   
	int string_length = 2;
	char c = 0;    
	FILE *file = 0;
	if(infile){
		if (!(file = fopen( infile, "r" ))){
			return 0;
 			fprintf(stderr,"Cannot open file '%s'\n", infile);
			exit(1);
		}
		if (fseek(file,0,SEEK_END) != 0){
			(void)fprintf(stderr, "ERROR: fseek failed\n");
			(void)exit(EXIT_FAILURE);
		}
		i= ftell (file);
		if (fseek(file,0,SEEK_START) != 0){
			(void)fprintf(stderr, "ERROR: fseek failed\n");
			(void)exit(EXIT_FAILURE);
		}
		string = malloc ((i+1)* sizeof(char));
		fread(string,sizeof(char), i, file);
		string[i] = 0;
		fclose(file);
	}else{
		if (!isatty(0)){
			string = malloc(sizeof(char*)*string_length);
			while (!feof (stdin)){
				c = getc(stdin);
				if (i == string_length){
					string_length <<=1;
					string = realloc(string,sizeof(char)*string_length);
				}
				string[i] = c;
				i++;
			}
			string[i-1] = 0;
		}else{
			return 0;
		}
	}
	return string;
}

struct alignment* read_sequences_from_swissprot(struct alignment* aln,char* string)
{
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	int i,j,c,n;
	char* p = 0;
	p = string;
	/*numseq = byg_count("ID   ",p);
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	aln = (struct alignment *) malloc(sizeof(struct alignment));
	numprofiles = (numseq << 1) - 1;
	aln->ft = 0;
	aln->si = 0;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);	
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}

	for (i = numseq;i--;){
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/
	c = 0;
	while(aln->sl[c]){
		c++;
	}

	
	while ((i = byg_end("ID   ",p)) != -1){
		p+=i;
		j = byg_start(" ",p);
		aln->lsn[c] = j;
		aln->sn[c] = malloc(sizeof(char)*(j+1));
		for (i = 0;i < j;i++){
			aln->sn[c][i] = p[i];
		}
		aln->sn[c][j] = 0;
		p+= j;
		j = byg_end("SQ   ",p);
		p+= j;
		j = byg_end("\n",p);
		p+= j;
		j = byg_start("//",p);

		aln->s[c] = malloc(sizeof(int)*(j+1));
		aln->seq[c] = malloc(sizeof(char)*(j+1));
		n = 0;
		for (i = 0;i < j;i++){
			if(isalpha((int)p[i])){
				aln->s[c][n] = aacode[toupper(p[i])-65];
				aln->seq[c][n] = p[i];
				n++;
			}
		}
		aln->s[c][n] = 0;
		aln->seq[c][n] = 0;
		aln->sl[c] = n;
		c++;
	}
	free(string);
	return aln;
}


struct alignment* read_alignment_from_swissprot(struct alignment* aln,char* string)
{
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	int i,j,c,n;
	char* p = 0;
	p = string;
	/*numseq = byg_count("ID   ",p);
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	aln = (struct alignment *) malloc(sizeof(struct alignment));
	numprofiles = (numseq << 1) - 1;
	aln->ft = 0;
	aln->si = 0;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);	
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}

	for (i = numseq;i--;){
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/
	c = 0;
	while(aln->sl[c]){
		c++;
	}

	fprintf(stderr,"found sequence:\n");
	while ((i = byg_end("ID   ",p)) != -1){
		p+=i;
		j = byg_start(" ",p);
		aln->lsn[c] = j;
		aln->sn[c] = malloc(sizeof(char)*(j+1));
		for (i = 0;i < j;i++){
			aln->sn[c][i] = p[i];
		}
		aln->sn[c][j] = 0;
		p+= j;
		j = byg_end("SQ   ",p);
		p+= j;
		j = byg_end("\n",p);
		p+= j;
		j = byg_start("//",p);
		fprintf(stderr,"found sequence:\n");
		aln->s[c] = malloc(sizeof(int)*(j+1));
		aln->seq[c] = malloc(sizeof(char)*(j+1));
		n = 0;
		for (i = 0;i < j;i++){
			if((int)p[i] > 32){
				if(isalpha((int)p[i])){
					aln->s[c][n] = aacode[toupper(p[i])-65];
				}else{
					aln->s[c][n] = -1;
				}
				fprintf(stderr,"%c",p[i]);
				aln->seq[c][n] = p[i];
				n++;
			}
		}
		
		fprintf(stderr,"\n\n");
		aln->s[c][n] = 0;
		aln->seq[c][n] = 0;
		aln->sl[c] = n;
		c++;
	}
	free(string);
	return aln;
}

struct alignment* read_sequences_macsim_xml(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	char *p = 0;
	int max = 0;
	
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

	/*aln = (struct alignment*) malloc(sizeof(struct alignment));
	numseq = byg_count("<seq-name>",string);
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	
	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->ft =  malloc(sizeof(struct feature* ) * (numseq));
	aln->si  =  malloc(sizeof(struct sequence_information* ) * (numseq));
	
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}
	for(i =0;i < numseq;i++){
		aln->ft[i] = 0;
		aln->si[i] = 0;
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/
		
	p = string;
	
	if(byg_count("<g>",p)){
		while((i = byg_start("<g>",p))!=-1){
			p+=i;
			j = byg_end("<r>",p);
			for(i = 0; i< j;i++){
				p[i] = ' ';
			}
			i = byg_start("</r>",p);
			p+=i;
			
			j = byg_end("</g>",p);
			for(i = 0; i< j;i++){
				p[i] = ' ';
			}
			
		}
	}
	p = string;

	c = 0;
	while(aln->sl[c]){
		c++;
	}
	

	
	while((i = byg_end("<sequence",p))!=-1){	
		p+=i;// p1 is at start of entry;
		max = byg_end("</sequence>",p);
			
		i = byg_end("<seq-name>",p);
		if(i < max){
			p +=i; //p1 is at the end of the sequence name tag
			j = byg_start("</seq-name>",p);
		
			aln->lsn[c] = j;
			aln->sn[c] = malloc(sizeof(char)*(j+1));
			for (i = 0;i < j;i++){
				aln->sn[c][i] = p[i];
			}
			aln->sn[c][j] = 0;
			
		}
		i = byg_end("<ftable>",p);
		if(i < max){
			aln->ft[c] = read_ft(aln->ft[c],p);
		}
		i = byg_end("<seq-data>",p);
		if(i < max){
			p+= i;
			j = byg_start("</seq-data>",p);
			aln->s[c] = malloc(sizeof(int)*(j+1));
			aln->seq[c] = malloc(sizeof(char)*(j+1));
			n = 0;
			for (i = 0;i < j;i++){
				if(isalpha((int)p[i])){
					aln->s[c][n] = aacode[toupper(p[i])-65];
					aln->seq[c][n] = p[i];
					n++;
				}
			}
			aln->s[c][n] = 0;
			aln->seq[c][n] = 0;
			aln->sl[c] = n;
		}
		
		c++;
	}
	free(string);
	return aln;
}


struct alignment* read_alignment_macsim_xml(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	char *p = 0;
	int max = 0;
	
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

	/*aln = (struct alignment*) malloc(sizeof(struct alignment));
	numseq = byg_count("<seq-name>",string);
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	
	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->ft =  malloc(sizeof(struct feature* ) * (numseq));
	aln->si  =  malloc(sizeof(struct sequence_information* ) * (numseq));
	
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}
	for(i =0;i < numseq;i++){
		aln->ft[i] = 0;
		aln->si[i] = 0;
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/
		
	p = string;
	
	if(byg_count("<g>",p)){
		while((i = byg_start("<g>",p))!=-1){
			p+=i;
			j = byg_end("<r>",p);
			for(i = 0; i< j;i++){
				p[i] = ' ';
			}
			i = byg_start("</r>",p);
			p+=i;
			
			j = byg_end("</g>",p);
			for(i = 0; i< j;i++){
				p[i] = ' ';
			}
			
		}
	}
	p = string;

	c = 0;
	while(aln->sl[c]){
		c++;
	}
	

	
	while((i = byg_end("<sequence",p))!=-1){	
		p+=i;// p1 is at start of entry;
		max = byg_end("</sequence>",p);
			
		i = byg_end("<seq-name>",p);
		if(i < max){
			p +=i; //p1 is at the end of the sequence name tag
			j = byg_start("</seq-name>",p);
		
			aln->lsn[c] = j;
			aln->sn[c] = malloc(sizeof(char)*(j+1));
			for (i = 0;i < j;i++){
				aln->sn[c][i] = p[i];
			}
			aln->sn[c][j] = 0;
			
		}
		i = byg_end("<ftable>",p);
		if(i < max){
			aln->ft[c] = read_ft(aln->ft[c],p);
		}
		i = byg_end("<seq-data>",p);
		if(i < max){
			p+= i;
			j = byg_start("</seq-data>",p);
			aln->s[c] = malloc(sizeof(int)*(j+1));
			aln->seq[c] = malloc(sizeof(char)*(j+1));
			n = 0;
			for (i = 0;i < j;i++){
				if((int)p[i]>32){
					if(isalpha((int)p[i])){
						aln->s[c][n] = aacode[toupper(p[i])-65];
					}else{
						aln->s[c][n] = -1;
					}
					aln->seq[c][n] = p[i];
					n++;
				}
			}
			aln->s[c][n] = 0;
			aln->seq[c][n] = 0;
			aln->sl[c] = n;
		}
		
		c++;
	}
	free(string);
	return aln;
}


struct feature* read_ft(struct feature* ft,char* p)
{

	int i,j;
	struct feature *n = 0;
	struct feature *old_n = 0;
	char tmp[10];
	char* p1 = 0;
	p1 = p;
	while((j = byg_end("<fitem>",p1))!= -1){
		i = byg_end("</seq-info>",p1);
		
		if(j >i){
			break;
		}

		n = malloc(sizeof(struct feature));
		n->next = 0;
		n->color = -1;

		p1+=j;// p1 is at start of entry;
		i = byg_end("<ftype>",p1);
		p1 +=i; //p1 is at the end of the sequence name tag
		j = byg_start("</ftype>",p1);

		n->type = malloc(sizeof(char*)*(j+1));
		for (i = 0; i < j;i++){
			n->type[i] = p1[i];
		}
		n->type[j] = 0;
		
		i = byg_end("<fstart>",p1);
		p1+= i;
		j = byg_start("</fstart>",p1);
		
		for (i = 0; i < j;i++){
			tmp[i] = p1[i];
		}
		tmp[j] = 0;
		n->start = atoi(tmp);
		i = byg_end("<fstop>",p1);
		p1+= i;
		j = byg_start("</fstop>",p1);
		for (i = 0; i < j;i++){
			tmp[i] = p1[i];
		}
		tmp[j] = 0;
		n->end = atoi(tmp);

		i = byg_end("<fnote>",p1);
		p1+= i;
		j = byg_start("</fnote>",p1);
		n->note = malloc(sizeof(char*)*(j+1));
		for (i = 0; i < j;i++){
			n->note[i] = p1[i];
		}
		
		n->note[j] = 0;

		
		if((old_n = ft)!= 0){
			while(old_n->next!=0){
				old_n = old_n->next;
			}
			old_n->next = n;
		}else{
			ft = n;
		}
		n = 0;
	}
	return ft;
}

struct alignment* read_sequences_uniprot_xml(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	char *p1 = 0;

	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

	/*aln = (struct alignment *) malloc(sizeof(struct alignment));
	numseq = byg_count("<entry",string);
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	
	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->si = 0;
	aln->ft = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}
	for(i =0;i < numseq;i++){
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/

	p1 = string;

	
	c = 0;
	while(aln->sl[c]){
		c++;
	}
	
	while((i = byg_end("<entry",p1))!=-1){

		p1+=i;// p1 is at start of entry;
		i = byg_end("<name>",p1);
		p1 +=i; //p1 is at the end of the sequence name tag
		j = byg_start("</name>",p1);
		aln->lsn[c] = j;
		aln->sn[c] = malloc(sizeof(char)*(j+1));
		for (i = 0;i < j;i++){
			aln->sn[c][i] = p1[i];
		}
		aln->sn[c][j] = 0;
		
		while((i = byg_end("<sequence",p1))!= -1 ){
			i = byg_end("<sequence",p1);
			p1+= i;
			i = byg_end(">",p1);
			p1 +=i;
		}
		
		j = byg_start("</sequence>",p1);

		aln->s[c] = malloc(sizeof(int)*(j+1));
		aln->seq[c] = malloc(sizeof(char)*(j+1));
		n = 0;
		for (i = 0;i < j;i++){
			if(isalpha((int)p1[i])){
				aln->s[c][n] = aacode[toupper(p1[i])-65];
				aln->seq[c][n] = p1[i];
				n++;
			}
		}
		aln->s[c][n] = 0;
		aln->seq[c][n] = 0;
		aln->sl[c] = n;
		c++;
	}
	free(string);
	return aln;
}



struct alignment* read_alignment_uniprot_xml(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	char *p1 = 0;

	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

	/*aln = (struct alignment *) malloc(sizeof(struct alignment));
	numseq = byg_count("<entry",string);
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	
	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->si = 0;
	aln->ft = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}
	for(i =0;i < numseq;i++){
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/

	p1 = string;

	
	c = 0;
	while(aln->sl[c]){
		c++;
	}
	
	while((i = byg_end("<entry",p1))!=-1){
		p1+=i;// p1 is at start of entry;
		i = byg_end("<name>",p1);
		p1 +=i; //p1 is at the end of the sequence name tag
		j = byg_start("</name>",p1);
		aln->lsn[c] = j;
		aln->sn[c] = malloc(sizeof(char)*(j+1));
		for (i = 0;i < j;i++){
			aln->sn[c][i] = p1[i];
		}
		aln->sn[c][j] = 0;
		i = byg_end("<sequence",p1);
		p1+= i;
		i = byg_end(">",p1);
		p1 +=i;
		j = byg_start("</sequence>",p1);
		aln->s[c] = malloc(sizeof(int)*(j+1));
		aln->seq[c] = malloc(sizeof(char)*(j+1));
		n = 0;
		for (i = 0;i < j;i++){
			if((int)p1[i] > 32){
				if(isalpha((int)p1[i])){
					aln->s[c][n] = aacode[toupper(p1[i])-65];
				}else{
					aln->s[c][n] = -1;
				}
				aln->seq[c][n] = p1[i];
				n++;
			}
		}
		aln->s[c][n] = 0;
		aln->seq[c][n] = 0;
		aln->sl[c] = n;
		c++;
	}
	free(string);
	return aln;
}

struct alignment* read_sequences_stockholm(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	char *p1 = 0;

	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

	/*aln = (struct alignment*) malloc(sizeof(struct alignment));
	p1 = string;
	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		if (!(byg_start("//",p1))){
			break;
		}	
		j = byg_end("#",p1);
		if(j != 1){
			numseq++;
		}
	}

	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->ft = 0;
	aln->si = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}
	for(i =0;i < numseq;i++){
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/

	c = 0;
	while(aln->sl[c]){
		c++;
	}

	p1 = string;
	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		if (!(byg_start("//",p1))){
			break;
		}	
		j = byg_end("#",p1);
		if(j != 1){
			j = byg_start(" ",p1);
			aln->lsn[c] = j;
			aln->sn[c] = malloc(sizeof(char)*(j+1));
			for (i = 0;i < j;i++){
				aln->sn[c][i] = p1[i];
			}
			aln->sn[c][j] = 0;

			
			p1+=j;
			j = byg_start("\n",p1);

			aln->s[c] = malloc(sizeof(int)*(j+1));
			aln->seq[c] = malloc(sizeof(char)*(j+1));
			n = 0;
			for (i = 0;i < j;i++){
				if(isalpha((int)p1[i])){
					aln->s[c][n] = aacode[toupper(p1[i])-65];
					aln->seq[c][n] = p1[i];
					n++;
				}
			}
			aln->s[c][n] = 0;
			aln->seq[c][n] = 0;
			aln->sl[c] = n;
			c++;
		}
	}

	free(string);
	return aln;
}

struct alignment* read_alignment_stockholm(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	char *p1 = 0;

	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};

	/*aln = (struct alignment*) malloc(sizeof(struct alignment));
	p1 = string;
	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		if (!(byg_start("//",p1))){
			break;
		}	
		j = byg_end("#",p1);
		if(j != 1){
			numseq++;
		}
	}

	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->ft = 0;
	aln->si = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}
	for(i =0;i < numseq;i++){
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
	}*/

	c = 0;
	while(aln->sl[c]){
		c++;
	}

	p1 = string;
	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		if (!(byg_start("//",p1))){
			break;
		}	
		j = byg_end("#",p1);
		if(j != 1){
			j = byg_start(" ",p1);
			aln->lsn[c] = j;
			aln->sn[c] = malloc(sizeof(char)*(j+1));
			for (i = 0;i < j;i++){
				aln->sn[c][i] = p1[i];
			}
			aln->sn[c][j] = 0;

			
			p1+=j;
			j = byg_start("\n",p1);

			aln->s[c] = malloc(sizeof(int)*(j+1));
			aln->seq[c] = malloc(sizeof(char)*(j+1));
			n = 0;
			for (i = 0;i < j;i++){
				if((int)p1[i] > 32){
					if(isalpha((int)p1[i])){
						aln->s[c][n] = aacode[toupper(p1[i])-65];
					}else{
						aln->s[c][n] = -1;
					}
					aln->seq[c][n] = p1[i];
					n++;
				}
			}
			aln->s[c][n] = 0;
			aln->seq[c][n] = 0;
			aln->sl[c] = n;
			c++;
		}
	}

	free(string);
	return aln;
}


struct alignment* read_sequences_clustal(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int len = 0;
	int i = 0;
	int j = 0;
	int start = 0;
	char *p1 = 0;
	int local_numseq = 0;

	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};


	//aln = (struct alignment*) malloc(sizeof(struct alignment));
	p1 = string;

	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		j = byg_end(" ",p1);
		n = byg_end("\n",p1);
		if(n > 2 && n>j && j!= 1){
			if(c ==0){
				i = j;
				while(p1[i] != '\n'){
					if (!isspace((int)p1[i])){
						len++;
					}
					i++;
				}		
			}
			c++;
		}else{
			if (c){
				if(c > local_numseq){
					local_numseq = c;
				}
				c =0;
			}
		}
	}

	/*numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->ft = 0;
	aln->si = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);

	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}

	for(i =0;i < numseq;i++){
		aln->lsn[i] = 0;
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
		aln->sl[i] = 0;*/
	start = 0;
	while(aln->sl[start]){
		start++;
	}

	for(i =start;i < local_numseq+start;i++){
		aln->s[i] = malloc(sizeof(int)*(len+1));
		aln->seq[i] = malloc(sizeof(char)*(len+1));
	}

	p1 = string;
	c = start;
	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		j = byg_end(" ",p1);
		n = byg_end("\n",p1);
		if(n > 2 && n>j && j!= 1){
			if(aln->lsn[c] == 0){
				aln->lsn[c] = j;
				aln->sn[c] = malloc(sizeof(char)*(j+1));
				for (i = 0;i < j;i++){
					aln->sn[c][i] = p1[i];
				}
				aln->sn[c][j] = 0;
			}
			for (i = j;i < n;i++){
				if(isalpha((int)p1[i])){
					aln->s[c][aln->sl[c]] = aacode[toupper(p1[i])-65];
					aln->seq[c][aln->sl[c]] = p1[i];
					aln->sl[c]++;
				}		
			}		
			c++;
		}else{
 			if (c != start){
				//c =0;
				c = start;
			}	
		}
	}
	for (i = start; i < local_numseq+start;i++){
		aln->s[i][aln->sl[i]] = 0;
	}
	free(string);
	return aln;
}


struct alignment* read_alignment_clustal(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int len = 0;
	int i = 0;
	int j = 0;
	int start = 0;
	char *p1 = 0;
	int local_numseq = 0;

	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	//int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,-1,13,14,15,16,17,-1,18,19,20,21,22};


	//aln = (struct alignment*) malloc(sizeof(struct alignment));
	p1 = string;

	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		j = byg_end(" ",p1);
		n = byg_end("\n",p1);
		if(n > 2 && n>j && j!= 1){
			if(c ==0){
				i = j;
				while(p1[i] != '\n'){
					if ((int)p1[i] > 32){
						len++;
					}
					i++;
				}		
			}
			c++;
		}else{
			if (c){
				if(c > local_numseq){
					local_numseq = c;
				}
				c =0;
			}
		}
	}

	/*numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq ));
	aln->seq = malloc(sizeof(char*) * (numseq ));
	aln->ft = 0;
	aln->si = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);

	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}

	for(i =0;i < numseq;i++){
		aln->lsn[i] = 0;
		aln->sip[i] = malloc(sizeof(int)*1);
		aln->nsip[i] = 1;
		aln->sip[i][0] = i;
		aln->sl[i] = 0;*/
	start = 0;
	while(aln->sl[start]){
		start++;
	}

	for(i =start;i < local_numseq+start;i++){
		aln->s[i] = malloc(sizeof(int)*(len+1));
		aln->seq[i] = malloc(sizeof(char)*(len+1));
	}

	p1 = string;
	c = start;
	while((i = byg_end("\n",p1))!=-1){
		p1+=i;
		j = byg_end(" ",p1);
		n = byg_end("\n",p1);
		if(n > 2 && n>j && j!= 1){
			if(aln->lsn[c] == 0){
				aln->lsn[c] = j;
				aln->sn[c] = malloc(sizeof(char)*(j+1));
				for (i = 0;i < j;i++){
					aln->sn[c][i] = p1[i];
				}
				aln->sn[c][j] = 0;
			}
			for (i = j;i < n;i++){
				if((int)p1[i] > 32){
					if(isalpha((int)p1[i])){
						aln->s[c][aln->sl[c]] = aacode[toupper(p1[i])-65];
					}else{
						aln->s[c][aln->sl[c]] = -1;
					}
					aln->seq[c][aln->sl[c]] = p1[i];
					aln->sl[c]++;
				}		
			}		
			c++;
		}else{
 			if (c != start){
				//c =0;
				c = start;
			}	
		}
	}
	for (i = start; i < local_numseq+start;i++){
		aln->s[i][aln->sl[i]] = 0;
		aln->seq[i][aln->sl[i]] = 0;
	}
	free(string);
	return aln;
}

struct alignment* read_sequences(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	int stop = 0;
	int start = 0;
	int nbytes;
	int local_numseq = 0;				// O	12				//U17
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	nbytes = strlen(string);

	//aln = (struct alignment*) malloc(sizeof(struct alignment));
	for (i =0;i < nbytes;i++){
		if (string[i] == '>'&& stop == 0){
			stop = 1;
			local_numseq++;
		}
		if (string[i] == '\n'){
			stop = 0;
		}
	}
	/*
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq));
	aln->seq = malloc(sizeof(char*) * (numseq));
	aln->ft = 0;
	aln->si = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}*/
	start = 0;
	while(aln->sl[start]){
		start++;
	}
	j = start;

	for (i =0;i < nbytes;i++){
		if (string[i] == '>' && stop == 0){
			stop = 1;
			aln->sl[j] =c;
			j++;
			c = 0;
		}
		if (string[i] == '\n'){
			if(stop == 1){
				aln->lsn[j-1] = n;
				n = 0;
			}
			stop = 0;
		}
		if (stop == 1 && string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
			n++;
		}
		if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
			if (isalpha((int)string[i])){
				c++;
			}
		}
	}
	aln->sl[j] = c;

	for (i =1+start;i < local_numseq+1+start;i++){
		if(!aln->sl[i]){
			fprintf(stderr,"Sequence %d has a length of 0!!\n",i-1);
			exit(1);
		}
		aln->sl[i-1] = aln->sl[i];
	}
	aln->sl[start+local_numseq] = 0;

	//for (i = numseq;i--;){
	for (i = start; i < local_numseq+start;i++){
		aln->s[i] = malloc(sizeof(int)*(aln->sl[i]+1));
		aln->seq[i] = malloc(sizeof(char)*(aln->sl[i]+1));
		aln->sn[i] = malloc(sizeof(char)*(aln->lsn[i]+1));
		//aln->sip[i] = malloc(sizeof(int)*1);
		//aln->nsip[i] = 1;
		//aln->sip[i][0] = i;
	}

	stop = 0;
	j = start;
	for (i =0;i < nbytes;i++){
		if (string[i] == '>' && stop == 0 ){
			stop = 1;
			j++;
			c = 0;
		}
		if (string[i] == '\n'){
			if(stop == 1){
				n = 0;
			}
			stop = 0;
		}
		if (stop == 1 &&string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
			aln->sn[j-1][n] = string[i];
			n++;
		}
		if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
			if(isalpha((int)string[i])){
				aln->s[j-1][c] = aacode[toupper(string[i])-65];
				aln->seq[j-1][c] = string[i];
				c++;
			}
		}
	}

	for (i = start;i< local_numseq+start;i++){
		aln->s[i][aln->sl[i]] = 0;
		aln->seq[i][aln->sl[i]] = 0;
		aln->sn[i][aln->lsn[i]] = 0;
	}	

	free(string);
	return aln;
}


struct alignment* read_alignment(struct alignment* aln,char* string)
{
	int c = 0;
	int n = 0;
	int i = 0;
	int j = 0;
	int stop = 0;
	int start = 0;
	int nbytes;
	int local_numseq = 0;				// O	12				//U17
	int aacode[26] = {0,1,2,3,4,5,6,7,8,-1,9,10,11,12,23,13,14,15,16,17,17,18,19,20,21,22};
	nbytes = strlen(string);

	//aln = (struct alignment*) malloc(sizeof(struct alignment));
	for (i =0;i < nbytes;i++){
		if (string[i] == '>'&& stop == 0){
			stop = 1;
			local_numseq++;
		}
		if (string[i] == '\n'){
			stop = 0;
		}
	}
	/*
	if(!numseq){
		fprintf(stderr,"No sequences found!\n");
		exit(1);
	}
	numprofiles = (numseq << 1) - 1;
	aln->s = malloc(sizeof(int*) * (numseq));
	aln->seq = malloc(sizeof(char*) * (numseq));
	aln->ft = 0;
	aln->si = 0;
	aln->sl = malloc(sizeof(int) * (numprofiles));
	aln->sip = malloc(sizeof(int*)* numprofiles);
	aln->nsip = malloc(sizeof(int)* numprofiles);
	aln->sn = malloc(sizeof(char*) * numseq);
	aln->lsn = malloc(sizeof(int) * numseq);
	
	for (i =0;i < numprofiles;i++){
		aln->sip[i] = 0;
		aln->nsip[i] = 0;
	}*/
	start = 0;
	while(aln->sl[start]){
		start++;
	}
	j = start;

	for (i =0;i < nbytes;i++){
		if (string[i] == '>' && stop == 0){
			stop = 1;
			aln->sl[j] =c;
			j++;
			c = 0;
		}
		if (string[i] == '\n'){
			if(stop == 1){
				aln->lsn[j-1] = n;
				n = 0;
			}
			stop = 0;
		}
		if (stop == 1 && string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
			n++;
		}
		if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
			if ((int)string[i] > 32){
				c++;
			}
		}
	}
	aln->sl[j] = c;

	for (i =1+start;i < local_numseq+1+start;i++){
		if(!aln->sl[i]){
			fprintf(stderr,"Sequence %d has a length of 0!!\n",i-1);
			exit(1);
		}
		aln->sl[i-1] = aln->sl[i];
	}
	aln->sl[start+local_numseq] = 0;
	//fprintf(stderr,"set to 0 : %d\n",start+local_numseq);
	//for (i = numseq;i--;){
	for (i = start; i < local_numseq+start;i++){
	//	fprintf(stderr,"len:%d %d\n",i,aln->sl[i]);
		aln->s[i] = malloc(sizeof(int)*(aln->sl[i]+1));
		aln->seq[i] = malloc(sizeof(char)*(aln->sl[i]+1));
		aln->sn[i] = malloc(sizeof(char)*(aln->lsn[i]+1));
		//aln->sip[i] = malloc(sizeof(int)*1);
		//aln->nsip[i] = 1;
		//aln->sip[i][0] = i;
	}

	stop = 0;
	j = start;
	for (i =0;i < nbytes;i++){
		if (string[i] == '>' && stop == 0 ){
			stop = 1;
			j++;
			c = 0;
		}
		if (string[i] == '\n'){
			if(stop == 1){
				n = 0;
			}
			stop = 0;
		}
		if (stop == 1 &&string[i] != '\n' && string[i] != 0 && string[i] != '>' ){
			aln->sn[j-1][n] = string[i];
			n++;
		}
		if (stop == 0 && string[i] != '\n' && string[i] != 0 ){
			if((int) string[i] > 32 ){
				if(isalpha((int)string[i])){
					aln->s[j-1][c] = aacode[toupper(string[i])-65];
				}else{
					aln->s[j-1][c] = -1;
				}
				aln->seq[j-1][c] = string[i];
				c++;
			}
		}
	}

	for (i = start;i< local_numseq+start;i++){
		aln->s[i][aln->sl[i]] = 0;
		aln->seq[i][aln->sl[i]] = 0;
		aln->sn[i][aln->lsn[i]] = 0;
	}	

	free(string);
	return aln;
}

