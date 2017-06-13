import descriptor

def feature( Sequence ):
	#!/usr/bin/python
	
	# Import module support
	#import support

	
	#Sequence = 'GPLGSGSKIKLEIYNETDMASASGYTPVPSVSEFQYIETETISNTPSPDLTVMSIDKSVLSPGESATITTIVKDIDGNPVNEVHINKTVARENLKGLWDYGPLKKENVPGKYTQVITYRGHSNERIDISFKYAMSFTKEISIRGRL'
	
	# Now you can call defined function that module as follows
	#support.print_func(Sequence)
	
	N = 6; # N=[0,5], if N=5 then use all descriptors
	   
	V=[];
	for j in range(0,N):
	    if j==1:
	        # Hydrophobicity
	        G1='RKEDQN'
	        G2='GASTPHY'
	        G3='CVLIMFW'
	        v = descriptor.descriptors(Sequence,G1,G2,G3);
	        for k in range(0,len(v)):
	            V.append(v[k])
	    if j==2:
	        # Normalized van der Waals volumns
	        G1=['G','A','S','C','T','P','D'];
	        G2=['N','V','E','Q','I','L'];
	        G3=['M','H','K','F','R','Y','W'];
	        v = descriptor.descriptors(Sequence,G1,G2,G3);
	        for k in range(0,len(v)):
	            V.append(v[k])
	    if j==3:
	        # Polarity
	        G1=['L','I','F','W','C','M','V','Y'];
	        G2=['P','A','T','G','S'];
	        G3=['H','Q','R','K','N','E','D'];
	        v = descriptor.descriptors(Sequence,G1,G2,G3);
	        for k in range(0,len(v)):
	            V.append(v[k])
	    if j==4:
	        # Polarizability
	        G1=['G','A','S','D','T'];
	        G2=['C','P','N','V','E','Q','A','L'];
	        G3=['K','M','H','F','R','Y','W'];
	        v = descriptor.descriptors(Sequence,G1,G2,G3);
	        for k in range(0,len(v)):
	            V.append(v[k])  
	    if j==5:
	        # Normalized frequency of alpha-helix
	        G1=['G','P','N','Y','C','S','T'];
	        G2=['R','H','D','V','W','I'];
	        G3=['Q','F','K','L','A','M','E'];
	        v = descriptor.descriptors(Sequence,G1,G2,G3);
	        for k in range(0,len(v)):
	            V.append(v[k])
	
	return V