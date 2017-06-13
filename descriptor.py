def descriptors( Sequence, G1, G2, G3 ):
    #print G1
    import math
    #print "Hello : ", G1, G2, G3
    Amino = 'ARNDCQEGHILKMFPSTWYV'
    
    S = Sequence;
    seq = []
    
    for j in range(0,len(S)):
        L1 = []
        L2 = []
        L3 = []
        for k in range(0,len(G1)):
            if G1[k]==S[j]:
                L1.append(1)
            else: L1.append(0)
        for k in range(0,len(G2)):
            if G2[k]==S[j]:
                L2.append(1)
            else: L2.append(0)
        for k in range(0,len(G3)):
            if G3[k]==S[j]:
                L3.append(1)
            else: L3.append(0)
        if sum(L1)==1:
            seq.append(1.0)
        elif sum(L2)==1:
            seq.append(2.0)
        elif sum(L3)==1:
            seq.append(3.0)
    
    
    length = len(seq);
    G = 3 # G is the number of groups
    
    # Occurence
    Occ = [];
    for j in range(0,len(Amino)):
        count = []
        for k in range(0,len(S)):
            if S[k]==Amino[j]:
                count.append(1.0)
            else: count.append(0.0)
        #print count
        Occ.append(sum(count)/length)
    #print Occ
    
    # Composition
    Comp = []
    for j in range(0,G):
        count = []
        for k in range(0,len(seq)):
            if seq[k]==(j+1):
                count.append(1.0)
            else: count.append(0.0)
        Comp.append(sum(count)/length)
    #print Comp
    
    # Transition
    t1 = []
    for k in range(0,len(seq)):
        if seq[k]==1:
            t1.append(1.0)
        else: t1.append(0.0)
    t2= []
    for k in range(0,len(seq)):
        if seq[k]==2:
            t2.append(1.0)
        else: t2.append(0.0)
    t3 = []
    for k in range(0,len(seq)):
        if seq[k]==3:
            t3.append(1.0)
        else: t3.append(0.0)

    #print t2
    #print t3
    T12=0; T13=0; T23=0;
    for j in range(0,length-1):
        s1=t1[j]+t1[j+1]
        s2=t2[j]+t2[j+1] #sum(t2[j:j+1]);
        s3=t3[j]+t3[j+1] #sum(t3[j:j+1]);
        #print s1, s2, s3
        if s1>=1 and s2>=1:
            T12 = T12 +1
        if s1>=1 and s3>=1:
            T13 = T13 +1
        if s2>=1 and s3>=1:
            T23 = T23 +1
    #print T12, T13, T23
    
    # Distribution
    P=[0.25,0.5,0.75,1];
    lenP=len(P);
    D1 = []
    D2 = []
    D3 = []
    

    c1= []
    r1= t1
    for m in range(0,len(t1)):
        c1.append(m)
    #Sorting descending order , with index
    for m in range(0,len(r1)-1):
        for n in range(m+1,len(r1)):
            if r1[n]>r1[m]:
                #print m, n
                temp = r1[n]
                r1[n]= r1[m]
                r1[m]= temp
                ct = c1[n]
                c1[n]= c1[m]
                c1[m]= ct
            
    c2= []
    r2= t2
    for m in range(0,len(t2)):
        c2.append(m)
    #Sorting descending order , with index
    for m in range(0,len(r2)-1):
        for n in range(m+1,len(r2)):
            if r2[n]>r2[m]:
                #print m, n
                temp = r2[n]
                r2[n]= r2[m]
                r2[m]= temp
                ct = c2[n]
                c2[n]= c2[m]
                c2[m]= ct
                
    c3= []
    r3= t3
    for m in range(0,len(t3)):
        c3.append(m)
    #Sorting descending order , with index
    for m in range(0,len(r3)-1):
        for n in range(m+1,len(r3)):
            if r3[n]>r3[m]:
                #print m, n
                temp = r3[n]
                r3[n]= r3[m]
                r3[m]= temp
                ct = c3[n]
                c3[n]= c3[m]
                c3[m]= ct
                
    if r1[1]!=0.0:
        D1.append(float(c1[1])/float(length))
    else: D1.append(0)
    if r2[1]!=0.0 :
        D2.append(float(c2[1])/float(length))
    else : D2.append(0)
    if r3[1]!=0.0 :
        D3.append(float(c3[1])/float(length))
    else: D3.append(0)
    
    
    
    for k in range(0,lenP):
        P1=int(math.ceil(sum(t1)*P[k]))
        if P1 != 0 :
            D1.append(float(c1[P1])/float(length))
        else:   D1[k+1]=0
        P2=int(math.ceil(sum(t2)*P[k]))
        if P2 != 0:
            D2.append(float(c2[P2])/float(length))
        else:   D2[k+1] = 0
        P3=int(math.ceil(sum(t3)*P[k]))
        if P3 != 0 :
            D3.append(float(c3[P3])/float(length))
        else: D3[k+1] = 0
    
    D = [D1,D2,D3]
  

    feature = []
    feature.append(Comp[0]); feature.append(Comp[2]); feature.append(Comp[2])
    feature.append(T12)
    feature.append(T13)
    feature.append(T23)
    feature.append(D1[0]); feature.append(D1[1]); feature.append(D1[2]); feature.append(D1[3]); feature.append(D1[4])
    feature.append(D2[0]); feature.append(D2[1]); feature.append(D2[2]); feature.append(D2[3]); feature.append(D2[4])
    feature.append(D3[0]); feature.append(D3[1]); feature.append(D3[2]); feature.append(D3[3]); feature.append(D3[4])
    
    
    return   feature