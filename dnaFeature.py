def extraction(pssmfile, spd3file):
    # print "Hello : ", par

    # !/usr/bin/python
    import cgi, cgitb
    cgitb.enable()

    import math
    import dubchak

    # import pydna
    # from pydna import extract_features

    # Create instance of FieldStorage
    form = cgi.FieldStorage()

    # Open PSSM file
    fo = open(pssmfile, "r+")
    # if fo.closed :
    #	fo = open("2LV4A.txt.pssm", "r+")
    str = fo.name + fo.read();
    # Close opend file
    fo.close()
    # str = "Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts A R N D C Q E G H I L K M F P S T W Y V A R N D C Q E G H I L K M F P S T W Y V 1 G 0 -2 0 -1 -2 -2 -2 5 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00 0.00 2 P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00 0.00 3 L -1 -2 -3 -4 -1 -2 -3 -4 -3 1 4 -2 2 0 -3 -2 -1 -2 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00 0.00 4 G 0 -2 0 -1 -2 -2 -2 5 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.00 0.00 5 "

    str = str.split()
    p = str[0:22]
    lastpos = str.index('Lambda')
    # print lastpos
    lastpos = lastpos - (lastpos % 62) - 4
    currentpos = str.index('Last') + 62
    p_seq = ''
    plen = 0
    pssm = {}
    while (currentpos < lastpos):
        p_no = [int(i) for i in str[currentpos]]
        p_seq = p_seq + str[currentpos + 1]
        # print "sequence = ",p_seq
        # print p_no[0], p_seq[p_no[0]-1]

        pssm[plen] = [int(i) for i in str[currentpos + 2:currentpos + 22]]
        currentpos = currentpos + 44
        plen = plen + 1
    # end while



    # Open SPD3 file
    fo = open(spd3file, "r+")
    # if fo.closed :
    #	fo = open("2LV4A.txt.spd3", "r+")
    spd_str = fo.name + fo.read();
    # Close opend file
    fo.close()
    spd_str = spd_str.split()
    currentpos = 10
    row = 0
    ss_seq = ''
    ASA = []
    phi = []
    psi = []
    theta = []
    tau = []
    coil = []
    betaSheet = []
    alphaHelix = []
    while (row < plen):
        ss_seq = ss_seq + spd_str[currentpos + 3]
        data = [float(i) for i in spd_str[currentpos + 4:currentpos + 12]]
        ASA.append(data[0])
        phi.append(data[1])
        psi.append(data[2])
        theta.append(data[3])
        tau.append(data[4])
        coil.append(data[5])
        betaSheet.append(data[6])
        alphaHelix.append(data[7])
        row = row + 1;
        currentpos = currentpos + 11
    # end while spd3


    # =============== Feature extraction start here======================
    # def extract_features(p_seq, plen, pssm):
    # feature 1: amino acid composition 20 features
    features = []
    F1 = []
    for i in range(0, 20):
        F1.append(0.0)
    AminoAcids = 'ARNDCQEGHILKMFPSTWYV'
    sequence = p_seq

    for i in range(0, plen):
        for j in range(0, 20):
            if sequence[i] == AminoAcids[j]:
                F1[j] = F1[j] + 1.0;  # Counting the individual peptide/aa occurance

    for j in range(0, 20):
        F1[j] = F1[j] / plen  # normalizing the occurance of aa with the seq-length
        features.append(F1[j])
    # return F1

    # F1 = extract_features(p_seq, plen, pssm)

    # ------------------------------------------
    F2 = dubchak.feature(sequence)
    for k in range(0, len(F2)):
        features.append(F2[k])
    # ------------------------------------------
    allmin = 999
    allmax = -999
    i = 0
    while (i < plen):
        rowmin = min(pssm[i][0:20])
        if rowmin < allmin:
            allmin = rowmin
        # print "rowmin=",rowmin
        rowmax = max(pssm[i][0:20])
        if rowmax > allmax:
            allmax = rowmax
        # print "rowmax=",rowmax
        i = i + 1
    # end while
    gmax = allmax
    gmin = allmin
    # print "gmax=",gmax
    # print "gmin=",gmin
    npssm = []
    # p = {}
    diff = float(gmax - gmin)
    # print "gmax-gmin=",diff
    for i in range(0, plen):
        a = []
        for j in range(0, 20):
            val = float(pssm[i][j] - gmin)
            # print val
            val = val / diff
            a.append(val)
        # print a
        npssm.append(a)
    # print "nppssm",npssm


    # feature 3: Bigram
    col = 20
    # plen = 102
    # bi-gram
    bigram_t = []
    bigram = []
    for i in range(0, col):
        for j in range(0, plen - 1):
            a = []
            if j + 1 < len(npssm):
                for k in range(0, col):
                    val = npssm[j][i] * npssm[j + 1][k]
                    a.append(val)
            bigram_t.append(a)
        bb = []
        for j in range(0, col):
            b = []
            for k in range(0, len(bigram_t) - 1):
                if k < len(bigram_t):
                    b.append(bigram_t[k][j])
            # print b
            # print "sum=",(sum(b)/plen)
            bb.append((sum(b) / plen))
        # bb = []
        bigram.append(bb)
    # print bb
    # print bigram
    # F3 = transpose(bigram)
    maxCol = len(bigram[0])
    for row in bigram:
        rowLength = len(row)
        if rowLength > maxCol:
            maxCol = rowLength
    F3 = []
    for colIndex in range(maxCol):
        F3.append([])
        for row in bigram:
            if colIndex < len(row):
                F3[colIndex].append(row[colIndex])
                features.append(row[colIndex])
    # print F3
    # -------------------------------------------------
    # feature 4: PSSM Composition
    pssm_occ = []
    for i in range(0, col):
        b = []
        for k in range(0, len(npssm)):
            if k < len(npssm):
                b.append(npssm[k][j])
                # print k
        pssm_occ.append((sum(b) / plen))
        features.append(sum(b) / plen)
    F4 = pssm_occ
    # print F4
    # -------------------------------------------------

    # feature 5:  PSSM Auto Covariance
    DF = 10
    # parameter for ACV
    acc_t = []
    acc = []
    for i in range(0, DF):
        for j in range(0, plen - i):
            a = []
            if j + 1 < len(npssm):
                for k in range(0, col):
                    val = npssm[j][k] * npssm[j + i][k]
                    a.append(val)
            acc_t.append(a)
        bb = []
        for j in range(0, col):
            b = []
            for k in range(0, len(acc_t) - 1):
                if k < len(acc_t):
                    b.append(acc_t[k][j])
            bb.append((sum(b) / plen))
        acc.append(bb)
        acc_t = []
    # print acc

    # Finding transpose of acc
    maxCol = len(acc[0])
    for row in acc:
        rowLength = len(row)
        if rowLength > maxCol:
            maxCol = rowLength
    F5 = []
    for colIndex in range(maxCol):
        F5.append([])
        for row in acc:
            if colIndex < len(row):
                F5[colIndex].append(row[colIndex])
                features.append(row[colIndex])
    # print F5

    # -------------------------------------------------
    # feature 6: 1-lead bigram
    one_l_bigram_t = []
    one_l_bigram = []
    for i in range(0, col):
        for j in range(0, plen - 2):
            a = []
            if j + 2 < len(npssm):
                for k in range(0, col):
                    val = npssm[j][i] * npssm[j + 2][k]
                    a.append(val)
            one_l_bigram_t.append(a)
        bb = []
        for j in range(0, col):
            b = []
            for k in range(0, len(one_l_bigram_t) - 1):
                if k < len(one_l_bigram_t):
                    b.append(one_l_bigram_t[k][j])
            bb.append((sum(b) / plen))
        one_l_bigram.append(bb)
        one_l_bigram_t = []

    maxCol = len(one_l_bigram[0])
    for row in one_l_bigram:
        rowLength = len(row)
        if rowLength > maxCol:
            maxCol = rowLength
    F6 = []
    for colIndex in range(maxCol):
        F6.append([])
        for row in one_l_bigram:
            if colIndex < len(row):
                F6[colIndex].append(row[colIndex])
                features.append(row[colIndex])

    # print F6
    # -------------------------------------------------
    # feature 7: segmented distribution
    # here we need a parameter Fp =5,10,25
    Fp = 10
    ACC = []
    k = 0;
    for j in range(0, 20):
        Tj = 0
        for m in range(0, len(npssm)):
            Tj = Tj + npssm[m][j]

        partialsum = 0;
        i = 0;
        tp = Fp
        while (tp <= 50):  # in range(Fp,50):
            tpTj = tp * Tj / 100
            while (partialsum <= tpTj and i <= len(npssm)):
                partialsum = partialsum + npssm[i][j]
                i = i + 1;
            ACC.append(i);
            features.append(i)
            k = k + 1;
            tp = tp + Fp
        # print k, tp
        partialsum = 0;
        i = len(npssm) - 1;
        index = 0;
        tp = Fp
        while (tp <= 50):  # for tp in range(Fp,50):
            while (partialsum <= tp * Tj / 100 and i >= 0):
                partialsum = partialsum + npssm[i][j];
                i = i - 1;
                index = index + 1;
            ACC.append(index);
            features.append(index)
            k = k + 1;
            tp = tp + Fp

    # print ACC
    F7 = ACC
    offf = k
    # -------------------------------------------------
    # feature 8: Secondary Structures composition Coil, Helix and Sheet
    F8 = []
    SS = 'CHE';
    for j in range(0, len(SS)):
        count = 0
        for i in range(0, len(ss_seq)):
            if ss_seq[i] == SS[j]:
                count = count + 1
        F8.append(count)
        features.append(count)
    # print F8

    # -------------------------------------------------
    F9 = []
    F9.append(float(F8[0]) / float(len(ss_seq)))
    F9.append(float(F8[1]) / float(len(ss_seq)))
    F9.append(float(F8[2]) / float(len(ss_seq)))

    features.append(float(F8[0]) / float(len(ss_seq)))
    features.append(float(F8[1]) / float(len(ss_seq)))
    features.append(float(F8[2]) / float(len(ss_seq)))
    # -------------------------------------------------
    # feature 10: ASA, angles occurence, probabilities CHE
    F10 = []
    pi = math.pi  # 3.141592653589793
    a = b = c = d = e = f = g = h = m = n = p = q = 0.0
    for i in range(0, len(ss_seq)):
        a = a + ASA[i];
        b = b + math.sin(phi[i] * pi / 180);
        c = c + math.cos(phi[i] * pi / 180);
        d = d + math.sin(psi[i] * pi / 180);
        e = e + math.cos(psi[i] * pi / 180);
        f = f + math.sin(theta[i] * pi / 180);
        g = g + math.cos(theta[i] * pi / 180);
        h = h + math.sin(tau[i] * pi / 180);
        m = m + math.cos(tau[i] * pi / 180);
        n = n + coil[i]
        p = p + betaSheet[i]
        q = q + alphaHelix[i]

    F10.append(a / len(ss_seq))
    F10.append(b / len(ss_seq))
    F10.append(c / len(ss_seq))
    F10.append(d / len(ss_seq))
    F10.append(e / len(ss_seq))
    F10.append(f / len(ss_seq))
    F10.append(g / len(ss_seq))
    F10.append(h / len(ss_seq))
    F10.append(m / len(ss_seq))
    F10.append(n / len(ss_seq))
    F10.append(p / len(ss_seq))
    F10.append(q / len(ss_seq))

    features.append(a / len(ss_seq))
    features.append(b / len(ss_seq))
    features.append(c / len(ss_seq))
    features.append(d / len(ss_seq))
    features.append(e / len(ss_seq))
    features.append(f / len(ss_seq))
    features.append(g / len(ss_seq))
    features.append(h / len(ss_seq))
    features.append(m / len(ss_seq))
    features.append(n / len(ss_seq))
    features.append(p / len(ss_seq))
    features.append(q / len(ss_seq))
    # -------------------------------------------------
    # feature 11: bigram of angles sine and cosine
    # first prepare a Lx8 array with sine cosine
    a_bigram_t = []
    a_bigram = []
    angles = []
    pi = math.pi
    a = []
    for i in range(0, len(ss_seq)):
        a.append(math.sin(phi[i] * pi / 180))
        a.append(math.cos(phi[i] * pi / 180))
        a.append(math.sin(psi[i] * pi / 180))
        a.append(math.cos(psi[i] * pi / 180))
        a.append(math.sin(theta[i] * pi / 180))
        a.append(math.cos(theta[i] * pi / 180))
        a.append(math.sin(tau[i] * pi / 180))
        a.append(math.cos(tau[i] * pi / 180))
        angles.append(a)
        a = []
    # print angles
    a = []
    for i in range(0, 8):
        for j in range(0, len(ss_seq) - 1):
            for k in range(0, 8):
                a.append(angles[j][i] * angles[j + 1][k])
                a_bigram_t.append(a)

        bb = []
        for j in range(0, 8):
            b = []
            for k in range(0, len(a_bigram_t) - 1):
                b.append(a_bigram_t[k][j])
            bb.append((sum(b) / len(ss_seq)))
        a_bigram.append(bb)

    # print a_bigram

    # Finding transpose of a_bigram
    maxCol = len(a_bigram[0])
    for row in a_bigram:
        rowLength = len(row)
        if rowLength > maxCol:
            maxCol = rowLength
    F11 = []
    for colIndex in range(maxCol):
        F11.append([])
        for row in a_bigram:
            if colIndex < len(row):
                F11[colIndex].append(row[colIndex])
                features.append(row[colIndex])
    # print F11
    # -------------------------------------------------
    # feature 12:  angles Auto Covariance
    DF = 10;  # parameter for ACV
    angles_acc_t = []
    angles_acc = []
    a = []
    for i in range(0, DF):
        for j in range(0, len(ss_seq) - i):
            for k in range(0, 8):
                a.append(angles[j][k] * angles[j + i][k])
            angles_acc_t.append(a)

        bb = []
        for j in range(0, 8):
            b = []
            for k in range(0, len(angles_acc_t) - 1):
                b.append(angles_acc_t[k][j])
            bb.append((sum(b) / len(ss_seq)))
        angles_acc.append(bb)
        angles_acc_t = [];

    # print angles_acc

    # Finding transpose of angles_acc
    maxCol = len(angles_acc[0])
    for row in angles_acc:
        rowLength = len(row)
        if rowLength > maxCol:
            maxCol = rowLength
    F12 = []
    for colIndex in range(maxCol):
        F12.append([])
        for row in angles_acc:
            if colIndex < len(row):
                F12[colIndex].append(row[colIndex])
                features.append(row[colIndex])
    # print F12
    # -------------------------------------------------
    # feature 13: bigram probabilities
    prob_bigram_t = []
    prob_bigram = []
    a = []
    for i in range(0, 3):
        for j in range(0, len(ss_seq) - 1):
            # for k in range(0,3):
            if i == 0:
                a.append(coil[j] * coil[j + 1])
                a.append(coil[j] * betaSheet[j + 1])
                a.append(coil[j] * alphaHelix[j + 1])
            if i == 1:
                a.append(betaSheet[j] * coil[j + 1])
                a.append(betaSheet[j] * betaSheet[j + 1])
                a.append(betaSheet[j] * alphaHelix[j + 1])
            if i == 2:
                a.append(alphaHelix[j] * coil[j + 1])
                a.append(alphaHelix[j] * betaSheet[j + 1])
                a.append(alphaHelix[j] * betaSheet[j + 1])
            prob_bigram_t.append(a)
            a = []

        bb = []
        for j in range(0, 3):
            b = []
            for k in range(0, len(prob_bigram_t) - 1):
                b.append(prob_bigram_t[k][j])
            bb.append((sum(b) / len(ss_seq)))
        prob_bigram.append(bb)
        prob_bigram_t = [];

    # print prob_bigram
    # Finding transpose of prob_bigram
    maxCol = len(prob_bigram[0])
    for row in prob_bigram:
        rowLength = len(row)
        if rowLength > maxCol:
            maxCol = rowLength
    F13 = []
    for colIndex in range(maxCol):
        F13.append([])
        for row in prob_bigram:
            if colIndex < len(row):
                F13[colIndex].append(row[colIndex])
                features.append(row[colIndex])
    # print F13

    # -------------------------------------------------
    # feature 14:  probabilities Auto Covariance
    DF = 10;  # parameter for ACV
    prob_acc_t = []
    prob_acc = []
    a = []
    for i in range(0, DF):
        for j in range(0, len(ss_seq) - i):
            a.append(coil[j] * coil[j + i])
            a.append(betaSheet[j] * betaSheet[j + i])
            a.append(alphaHelix[j] * alphaHelix[j + i])  # need to check again
            prob_acc_t.append(a)
            a = []

        bb = []
        for j in range(0, 3):
            b = []
            for k in range(0, len(prob_acc_t) - 1):
                b.append(prob_acc_t[k][j])
            bb.append((sum(b) / len(ss_seq)))
        prob_acc.append(bb)
        prob_acc_t = []

    # print prob_acc

    # Finding transpose of prob_acc
    maxCol = len(prob_acc[0])
    for row in prob_acc:
        rowLength = len(row)
        if rowLength > maxCol:
            maxCol = rowLength
    F14 = []
    for colIndex in range(maxCol):
        F14.append([])
        for row in prob_acc:
            if colIndex < len(row):
                F14[colIndex].append(row[colIndex])
                features.append(row[colIndex])

    # print F14
    # -----------------------------------------------------------------
    # Concat all features in one list


    # -----------------------------------------------------------------


    # =============== Feature extraction END here======================



    return features