thesis.pdf: thesis.dvi 
	dvipdfmx thesis.dvi

thesis.dvi: thesis.tex 1_intro.tex 2_clay.tex 3_theory.tex 4_experiment.tex\
	Figs/fig2_1.eps\
	Figs/fig3_1.eps Figs/fig3_2.eps Figs/fig3_3.eps Figs/fig3_4.eps Figs/fig3_5.eps Figs/fig3_6.eps Figs/fig3_7.eps Figs/fig3_8.eps\
	Figs/fig4_1.eps Figs/fig4_2.eps Figs/fig4_3.eps Figs/fig4_4.eps
	platex thesis.tex

Figs/fig2_1.eps: Figs/smectite.svgz
	inkscape -z -f Figs/smectite.svgz -E Figs/fig2_1.eps
Figs/fig3_1.eps: Figs/R.svgz
	inkscape -z -f Figs/R.svgz -E Figs/fig3_1.eps
Figs/fig3_2.eps: Figs/C.svgz
	inkscape -z -f Figs/C.svgz -E Figs/fig3_2.eps
Figs/fig3_3.eps: Figs/L.svgz
	inkscape -z -f Figs/L.svgz -E Figs/fig3_3.eps
Figs/fig3_4.eps: Figs/RC.svgz
	inkscape -z -f Figs/RC.svgz -E Figs/fig3_4.eps
Figs/fig3_5.eps: Figs/CC.svgz
	inkscape -z -f Figs/CC.svgz -E Figs/fig3_5.eps
Figs/fig3_6.eps: Figs/CD.svgz
	inkscape -z -f Figs/CD.svgz -E Figs/fig3_6.eps
Figs/fig3_7.eps: Figs/CPE.svgz
	inkscape -z -f Figs/CPE.svgz -E Figs/fig3_7.eps
Figs/fig3_8.eps: Figs/Wbrg.svgz
	inkscape -z -f Figs/Wbrg.svgz -E Figs/fig3_8.eps

Figs/fig4_1.eps: Figs/sample.svgz
	inkscape -z -f Figs/sample.svgz -E Figs/fig4_1.eps
Figs/fig4_2.eps: Figs/setup.svgz
	inkscape -z -f Figs/setup.svgz -E Figs/fig4_2.eps
Figs/fig4_3.eps: Figs/cabling.svgz
	inkscape -z -f Figs/cabling.svgz -E Figs/fig4_3.eps
Figs/fig4_4.eps: Figs/rhod.svgz
	inkscape -z -f Figs/rhod.svgz -E Figs/fig4_4.eps
