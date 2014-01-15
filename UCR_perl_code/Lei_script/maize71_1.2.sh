2bwt-builder miRNAInGnm/maize/maize.chr01.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr01.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr01.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr01.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr01.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr01.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr01.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr01.M0r2v0.soap miRNAInGnm/maize/maize.chr01.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr01.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr01.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr01.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr01.n1.out_71 2>_maize10.2cn20.CDS.M0.chr01.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr02.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr02.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr02.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr02.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr02.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr02.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr02.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr02.M0r2v0.soap miRNAInGnm/maize/maize.chr02.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr02.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr02.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr02.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr02.n1.out_71 2>_maize10.2cn20.CDS.M0.chr02.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr03.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr03.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr03.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr03.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr03.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr03.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr03.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr03.M0r2v0.soap miRNAInGnm/maize/maize.chr03.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr03.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr03.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr03.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr03.n1.out_71 2>_maize10.2cn20.CDS.M0.chr03.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr04.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr04.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr04.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr04.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr04.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr04.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr04.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr04.M0r2v0.soap miRNAInGnm/maize/maize.chr04.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr04.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr04.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr04.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr04.n1.out_71 2>_maize10.2cn20.CDS.M0.chr04.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr05.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr05.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr05.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr05.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr05.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr05.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr05.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr05.M0r2v0.soap miRNAInGnm/maize/maize.chr05.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr05.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr05.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr05.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr05.n1.out_71 2>_maize10.2cn20.CDS.M0.chr05.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr06.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr06.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr06.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr06.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr06.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr06.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr06.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr06.M0r2v0.soap miRNAInGnm/maize/maize.chr06.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr06.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr06.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr06.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr06.n1.out_71 2>_maize10.2cn20.CDS.M0.chr06.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr07.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr07.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr07.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr07.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr07.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr07.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr07.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr07.M0r2v0.soap miRNAInGnm/maize/maize.chr07.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr07.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr07.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr07.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr07.n1.out_71 2>_maize10.2cn20.CDS.M0.chr07.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr08.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr08.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr08.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr08.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr08.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr08.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr08.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr08.M0r2v0.soap miRNAInGnm/maize/maize.chr08.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr08.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr08.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr08.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr08.n1.out_71 2>_maize10.2cn20.CDS.M0.chr08.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr09.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr09.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr09.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr09.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr09.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr09.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr09.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr09.M0r2v0.soap miRNAInGnm/maize/maize.chr09.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr09.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr09.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr09.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr09.n1.out_71 2>_maize10.2cn20.CDS.M0.chr09.n1.err_71
2bwt-builder miRNAInGnm/maize/maize.chr10.fa
soap -a miRNAInGnm/maize_sRNAs.fa -D miRNAInGnm/maize/maize.chr10.fa.index -M 0 -r 2 -v 0 -o miRNAInGnm/maize/maize_sRNAs.VS.Zmachr10.M0r2v0.soap 2> miRNAInGnm/maize/maize_sRNAs.VS.Zmachr10.M0r2v0.soap.out
perl miRNAInGnm/soapFilter.pl -0 0 -i miRNAInGnm/maize/maize_sRNAs.VS.Zmachr10.M0r2v0.soap -p miRNAInGnm/maize/maize10-rep.VS.Chr10.r2v0.2cn20.CDS.soap > miRNAInGnm/maize/maize_sRNAs.VS.Zmachr10.M0r2v0.filter.soap
rm miRNAInGnm/maize/maize_sRNAs.VS.Zmachr10.M0r2v0.soap miRNAInGnm/maize/maize.chr10.fa.index*
perl miRNAInGnm/miRNAInGnm7.1s.pl -p miRNAInGnm/ -i miRNAInGnm/maize/maize10-rep.VS.Chr10.r2v0.2cn20.CDS.soap -f 0 -d miRNAInGnm/maize_sRNAs.fa -g miRNAInGnm/maize/maize.chr10.fa -S miRNAInGnm/maize/maize_sRNAs.VS.Zmachr10.M0r2v0.filter.soap -q "-b \"\"" -s 1 -c "-m 4 -n 1 -l 12100 -R 0.75 -a 5 -M 0" -m 45000 -e "-n,1,-r,0.75" -j soap -0 0 -a 0 >_maize10.2cn20.CDS.M0.chr10.n1.out_71 2>_maize10.2cn20.CDS.M0.chr10.n1.err_71

