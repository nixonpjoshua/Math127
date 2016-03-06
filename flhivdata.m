% flhivdata.m
% 
% Sequence data from Florida dentist AIDS cluster, from GenBank, is stored in 
% variables: dnt, lc1, lc5, ptb, ptc, ptd. A 20-site excerpt of all of these 
% is stored in the variable: parsimony.
%
% Note: Comments show full GenBank entries. Sequences are not of equal length, 
% but they are aligned from start. 
%
% 8/2/03


%LOCUS       HIVFLD1       360 bp ss-RNA             VRL       09-JUL-1992
%DEFINITION  Human immunodeficiency virus type 1, viral sample FLD1, V3 region.
%ACCESSION   M90848
%KEYWORDS    .
%SOURCE      Human immunodeficiency virus type 1 (HIV-1), M13 clone 1 of DNA ID
%            4430.
%  ORGANISM  Human immunodeficiency virus type 1
%            Viridae; ss-RNA enveloped viruses; Positive strand RNA virus;
%            Retroviridae; Lentivirinae.
%REFERENCE   1  (bases 1 to 360)
%  AUTHORS   Ou,C.-Y.-., Ciesielski,C.A., Myers,G., Bandea,C.I., Luo,C.C.,
%            Korber,B.T., Mullins,J.I., Schochetman,G., Berkelman,R.L.,
%            Economou,A.N., Witte,J.J., Furman,L.J., Satten,G.A., Curran,J.W.
%            and Jaffe,H.W.
%  TITLE     Molecular Epidemiology of HIV Transmission in a Dental Practice
%  JOURNAL   Science 256, 1165-1171 (1992)
%  STANDARD  full automatic
%COMMENT     Kindly submitted in computer readable form by the CDC (Centers for
%            Disease Control), Atlanta, GA.
%            The sequence in this entry is one of 6 clone sequences over the V3
%            region obtained by the CDC from a Florida dentist.
%            Please note that for this set of sequences, clone numbers from the
%            V3 region do not correspond with similar numbers from the V4C3V5
%            region.
%            NCBI gi: 326847
%FEATURES             Location/Qualifiers
%     CDS             <1..>360
%                     /note="env polyprotein;  NCBI gi: 326848."
%                     /codon_start=1
%                     /translation="LAEEEVVIRSANFTDNAKIIIVQLNASVEINCTRPNNNTRKGIH
%                     IGPGRAFYATGEIIGDIRQAHCNISREKWNNTLNQVVTELREQFGNKTITFNHSSGGD
%                     PEIVMHSFNCGGEFFYCN"
%     source          1..360
%                     /organism="Human immunodeficiency virus type 1"
%BASE COUNT      148 a     50 c     73 g     89 t
%ORIGIN      near the C2-V3 boundary
d0='ctagcagaagaagaggtagtaattagatctgccaatttcacagacaatgctaaaatcataatagtacagctgaatgcatctgtagaaattaattgtacaagg';
d1='cccaacaacaatacaagaaaaggtatacatataggaccagggagagcattttatgcaacaggagaaataataggagatataagacaagcacattgtaacatt';
d2='agtagagaaaaatggaataatactttaaaccaggtagttacagaattaagggaacaatttgggaataaaacaataacctttaatcactcctcaggaggggacccagaaattgtaatgcacagttttaattgtggaggggaatttttctattgtaat';
dnt=[d0 d1 d2];
clear d0 d1 d2;

%LOCUS       HIVFLPBD3     343 bp ss-RNA             VRL       09-JUL-1992
%DEFINITION  Human immunodeficiency virus type 1, viral sample FLPBD3, V3
%            region.
%ACCESSION   M90873
%KEYWORDS    .
%SOURCE      Human immunodeficiency virus type 1 (HIV-1), M13 clone DA3 of DNA
%            ID 5592.
%  ORGANISM  Human immunodeficiency virus type 1
%            Viridae; ss-RNA enveloped viruses; Positive strand RNA virus;
%            Retroviridae; Lentivirinae.
%REFERENCE   1  (bases 1 to 343)
%  AUTHORS   Ou,C.-Y.-., Ciesielski,C.A., Myers,G., Bandea,C.I., Luo,C.C.,
%            Korber,B.T., Mullins,J.I., Schochetman,G., Berkelman,R.L.,
%            Economou,A.N., Witte,J.J., Furman,L.J., Satten,G.A., Curran,J.W.
%            and Jaffe,H.W.
%  TITLE     Molecular Epidemiology of HIV Transmission in a Dental Practice
%  JOURNAL   Science 256, 1165-1171 (1992)
%  STANDARD  full automatic
%COMMENT     Kindly submitted in computer readable form by the CDC (Centers for
%            Disease Control), Atlanta, GA.
%            The sequence in this entry is one of 13 clone sequences over the V3
%            region obtained by the CDC from patient 'B' of a Florida dentist.
%            Please note that for this set of sequences, clone numbers from the
%            V3 region do not correspond with similar numbers from the V4C3V5
%            region.
%            NCBI gi: 326972
%FEATURES             Location/Qualifiers
%     CDS             <1..>343
%                     /note="env polyprotein;  NCBI gi: 326973."
%                     /codon_start=1
%                     /translation="LAEEEIVIRSANFTDNAKIIIVQLNASVEINCTRPDNNTRKGIH
%                     IGPGRAFYATGEIIGDIRQAHCNISGAKWNNTIEQVKTKLREQFGNKTIIFNHSSGGD
%                     PEIVMHSFNCGGX"
%     source          1..343
%                     /organism="Human immunodeficiency virus type 1"
%BASE COUNT      146 a     48 c     72 g     77 t
%ORIGIN      near the C2-V3 boundary
ptb1='ctagcagaagaagagatagtaattagatctgccaatttcacagacaatgctaaaatcata';
ptb2='atagtacagctgaatgcatctgtagaaattaattgtacaagacccgacaacaatacaaga';
ptb3='aaaggtatacatataggaccagggagggcattttatgcaacaggagaaataataggagat';
ptb4='ataagacaagcacattgtaacattagtggagcaaaatggaataatactatagaacaggta';
ptb5='aagacaaaattaagagaacaatttgggaataaaacaataatctttaatcactcctcagga';
ptb6='ggggacccagaaattgtaatgcacagttttaattgtggagggg';
ptb=[ptb1 ptb2 ptb3 ptb4 ptb5 ptb6];
clear ptb1 ptb2 ptb3 ptb4 ptb5 ptb6;

%LOCUS       HIVFLPC14     349 bp ss-RNA             VRL       09-JUL-1992
%DEFINITION  Human immunodeficiency virus type 1, viral sample FLPC14, V3
%            region.
%ACCESSION   M90878
%KEYWORDS    .
%SOURCE      Human immunodeficiency virus type 1 (HIV-1), M13 clone A14 of DNA
%            ID 5606.
%  ORGANISM  Human immunodeficiency virus type 1
%            Viridae; ss-RNA enveloped viruses; Positive strand RNA virus;
%            Retroviridae; Lentivirinae.
%REFERENCE   1  (bases 1 to 349)
%  AUTHORS   Ou,C.-Y.-., Ciesielski,C.A., Myers,G., Bandea,C.I., Luo,C.C.,
%            Korber,B.T., Mullins,J.I., Schochetman,G., Berkelman,R.L.,
%            Economou,A.N., Witte,J.J., Furman,L.J., Satten,G.A., Curran,J.W.
%            and Jaffe,H.W.
%  TITLE     Molecular Epidemiology of HIV Transmission in a Dental Practice
%  JOURNAL   Science 256, 1165-1171 (1992)
%  STANDARD  full automatic
%COMMENT     Kindly submitted in computer readable form by the CDC (Centers for
%            Disease Control), Atlanta, GA.
%            The sequence in this entry is one of 5 clones over the V3 region
%            obtained by the CDC from patient 'C' of a Florida dentist.
%            Please note that for this set of sequences, clone numbers from the
%            V3 region do not correspond with similar numbers from the V4C3V5
%            region.
%            NCBI gi: 327000
%FEATURES             Location/Qualifiers
%     CDS             <1..>349
%                     /note="env polyprotein;  NCBI gi: 327001."
%                     /codon_start=1
%                     /translation="LAEEEVVIRSADFTDNAKIIIVQLNASVEINCTRPNNNTRKGIH
%                     IGPGRAVYATDRIIGDIRQAHCNISREKWNNTLKQVVTKLREQFVNKTIIFTHPSGGD
%                     PEIVMHSVNCGGEFX"
%     source          1..349
%                     /organism="Human immunodeficiency virus type 1"
%BASE COUNT      148 a     49 c     69 g     83 t
%ORIGIN      near the C2-V3 boundary
ptc1='ctagcagaagaagaggtagtaattagatctgccgatttcacagacaatgctaaaatcataatagtacagctaaatgcatctgtagaaattaattgtacaagacctaacaacaatac';
ptc2='aagaaaaggtatacatataggaccagggagagcagtttatgcaacagacagaataataggagatataagacaagcacattgtaacattagtagagaaaaatggaataatactttaaaacagg';
ptc3='tagttacaaaattaagagaacaatttgtgaataaaacaataatctttactcacccctcaggaggggacccagaaattgtaatgcacagtgttaattgtggaggggaatttt';
ptc=[ptc1 ptc2 ptc3];
clear ptc1 ptc2 ptc3;

%//
%LOCUS       HIVFLPD9      337 bp ss-RNA             VRL       09-JUL-1992
%DEFINITION  Human immunodeficiency virus type 1, viral sample FLPD9, V3 region.
%ACCESSION   M90885
%KEYWORDS    .
%SOURCE      Human immunodeficiency virus type 1 (HIV-1), M13 clone A09 of DNA
%            ID 5826.
%  ORGANISM  Human immunodeficiency virus type 1
%            Viridae; ss-RNA enveloped viruses; Positive strand RNA virus;
%            Retroviridae; Lentivirinae.
%REFERENCE   1  (bases 1 to 337)
%  AUTHORS   Ou,C.-Y.-., Ciesielski,C.A., Myers,G., Bandea,C.I., Luo,C.C.,
%            Korber,B.T., Mullins,J.I., Schochetman,G., Berkelman,R.L.,
%            Economou,A.N., Witte,J.J., Furman,L.J., Satten,G.A., Curran,J.W.
%            and Jaffe,H.W.
%  TITLE     Molecular Epidemiology of HIV Transmission in a Dental Practice
%  JOURNAL   Science 256, 1165-1171 (1992)
%  STANDARD  full automatic
%COMMENT     Kindly submitted in computer readable form by the CDC (Centers for
%            Disease Control), Atlanta, GA.
%            The sequence in this entry is one of 5 clones over the V3 region
%            obtained by the CDC from patient 'D' of a Florida dentist.
%            Please note that for this set of sequences, clone numbers from the
%            V3 region do not correspond with similar numbers from the V4C3V5
%            region.
%            NCBI gi: 327038
%FEATURES             Location/Qualifiers
%     CDS             <1..>337
%                     /note="env polyprotein;  NCBI gi: 327039."
%                     /codon_start=1
%                     /translation="LAEEEVVIRSANFSDNAKTIIVQLNKSVKIPCIRPSNNTRQSIP
%                     IGPGKAVYATGQIIGDIRKAHRNLSEAIWNNTLKQIVKKLKEQFKNKTIVFNQSSGGD
%                     PEIVMHSFNCX"
%     source          1..337
%                     /organism="Human immunodeficiency virus type 1"
%BASE COUNT      146 a     51 c     62 g     78 t
%ORIGIN      near the C2-V3 boundary
ptd1='ctagcagaagaagaggtagtaattagatctgcaaatttctcggacaatgctaaaaccataatagtacagctgaataaatctgtaaaaattccttgtataagacccag';
ptd2='caataatacaagacaaagtatacctataggaccagggaaagcagtttatgcaacaggacagataataggagatataagaaaggcacatcgtaaccttagtgaagcaatatggaataacacgttaaaacagatagttaaaaaat';
ptd3='taaaagaacaatttaagaataaaacaatagtcttcaatcaatcctcaggaggggacccagaaattgtaatgcacagttttaattgtg';
ptd=[ptd1 ptd2 ptd3];
clear ptd1 ptd2 ptd3;

%LOCUS       HIVFLQ511     325 bp ss-RNA             VRL       09-JUL-1992
%DEFINITION  Human immunodeficiency virus type 1, viral sample LC01.DA11, V3
%            region.
%ACCESSION   M90915
%KEYWORDS    .
%SOURCE      Human immunodeficiency virus type 1 (HIV-1), M13 clone DA11 of DNA
%            ID 5240.
%  ORGANISM  Human immunodeficiency virus type 1
%            Viridae; ss-RNA enveloped viruses; Positive strand RNA virus;
%            Retroviridae; Lentivirinae.
%REFERENCE   1  (bases 1 to 325)
%  AUTHORS   Ou,C.-Y.-., Ciesielski,C.A., Myers,G., Bandea,C.I., Luo,C.C.,
%            Korber,B.T., Mullins,J.I., Schochetman,G., Berkelman,R.L.,
%            Economou,A.N., Witte,J.J., Furman,L.J., Satten,G.A., Curran,J.W.
%            and Jaffe,H.W.
%  TITLE     Molecular Epidemiology of HIV Transmission in a Dental Practice
%  JOURNAL   Science 256, 1165-1171 (1992)
%  STANDARD  full automatic
%COMMENT     Kindly submitted in computer readable form by the CDC (Centers for
%            Disease Control), Atlanta, GA.
%            The sequence in this entry is one of 2 clones over the V3 region
%            obtained by the CDC from this Florida control sample.
%            Please note that for this set of sequences, clone numbers from the
%            V3 region do not correspond with similar numbers from the V4C3V5
%            region.
%            NCBI gi: 327201
%FEATURES             Location/Qualifiers
%     CDS             <1..>325
%                     /note="env polyprotein;  NCBI gi: 327202."
%                     /codon_start=1
%                     /translation="LAEEEVVIRSENFTNNAKIIIVHLNKTVNITCTRPNNNTRRSIP
%                     MGPGKAFYTTEIIGNIRQAHCNLSKAEWNNTLRQIVKKLREQFKNKTIVFNHSSGGDP
%                     EIVMHSX"
%     source          1..325
%                     /organism="Human immunodeficiency virus type 1"
%BASE COUNT      148 a     49 c     54 g     74 t
%ORIGIN      near the C2-V3 boundary
lc1a='ctagcagaagaagaagtagtaattagatctgaaaatttcacgaataatgctaaaatcataatagtacacctgaataaaactgtaaatattacttgtacaagacccaacaacaatacaagaaga';
lc1b='agtatacctatgggaccagggaaagcattttatacaacagaaataataggaaatataagacaagcacattgtaaccttagtaaagcagaatggaataacactttaagacagatagttaaaaagttaa';
lc1c='gagaacaatttaagaataaaacaatagtcttcaatcactcctcaggaggggacccagaaattgtaatgcacagtt';
lc1=[lc1a lc1b lc1c];
clear lc1a lc1b lc1c;

%//
%LOCUS       HIVFLQ11      330 bp ss-RNA             VRL       09-JUL-1992
%DEFINITION  Human immunodeficiency virus type 1, viral sample LC05.EA14N, V3
%            region.
%ACCESSION   M90934
%KEYWORDS    .
%SOURCE      Human immunodeficiency virus type 1 (HIV-1), direct PCR
%            amplification product of Florida control sample EA14N of DNA ID
%            5292.
%  ORGANISM  Human immunodeficiency virus type 1
%            Viridae; ss-RNA enveloped viruses; Positive strand RNA virus;
%            Retroviridae; Lentivirinae.
%REFERENCE   1  (bases 1 to 330)
%  AUTHORS   Ou,C.-Y.-., Ciesielski,C.A., Myers,G., Bandea,C.I., Luo,C.C.,
%            Korber,B.T., Mullins,J.I., Schochetman,G., Berkelman,R.L.,
%            Economou,A.N., Witte,J.J., Furman,L.J., Satten,G.A., Curran,J.W.
%            and Jaffe,H.W.
%  TITLE     Molecular Epidemiology of HIV Transmission in a Dental Practice
%  JOURNAL   Science 256, 1165-1171 (1992)
%  STANDARD  full automatic
%COMMENT     Kindly submitted in computer readable form by the CDC (Centers for
%            Disease Control), Atlanta, GA.
%            Please note that for this set of sequences, clone numbers from the
%            V3 region do not correspond with similar numbers from the V4C3V5
%            region.
%            NCBI gi: 327131
%FEATURES             Location/Qualifiers
%     CDS             <1..>330
%                     /note="env polyprotein;  NCBI gi: 327132."
%                     /codon_start=1
%                     /translation="LAEEEVVIRSENFTNNAKTIIVQLKESVKINCIRPNNNTRRSIN
%                     MGPGRAFYTTGDIIGDIRQAHCNISKAEWNNTLKQIVQKLKEQFRNKTIVFNQSSGGD
%                     PEVVTHSF"
%     source          1..330
%                     /organism="Human immunodeficiency virus type 1"
%BASE COUNT      149 a     46 c     61 g     74 t
%ORIGIN      near the C2-V3 boundary
lc5a='ctagcagaagaggaggtagtaattagatctgaaaatttcacgaacaatgctaaaaccata';
lc5b='atagtacaactgaaagaatctgtaaaaattaattgtataagacccaacaacaatacaaga';
lc5c='agaagtataaatatgggaccagggagggcattttatacaacaggagacataataggagat';
lc5d='ataagacaagcacattgtaacattagtaaagcagagtggaataacactttaaaacagata';
lc5e='gttcaaaaattaaaagaacaatttaggaataaaacaatagtctttaatcaatcctcagga';
lc5f='ggggacccagaagttgtaacacacagtttt';
lc5=[lc5a lc5b lc5c lc5d lc5e lc5f];
clear lc5a lc5b lc5c lc5d lc5e lc5f;

%convert to upper case

dnt(dnt=='a')='A';dnt(dnt=='c')='C';dnt(dnt=='g')='G';dnt(dnt=='t')='T';
lc1(lc1=='a')='A';lc1(lc1=='c')='C';lc1(lc1=='g')='G';lc1(lc1=='t')='T';
lc5(lc5=='a')='A';lc5(lc5=='c')='C';lc5(lc5=='g')='G';lc5(lc5=='t')='T';
ptb(ptb=='a')='A';ptb(ptb=='c')='C';ptb(ptb=='g')='G';ptb(ptb=='t')='T';
ptc(ptc=='a')='A';ptc(ptc=='c')='C';ptc(ptc=='g')='G';ptc(ptc=='t')='T';
ptd(ptd=='a')='A';ptd(ptd=='c')='C';ptd(ptd=='g')='G';ptd(ptd=='t')='T';





