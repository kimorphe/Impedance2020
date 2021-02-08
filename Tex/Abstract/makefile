main.pdf: main.dvi
	dvipdfmx -p a4 main.dvi
main.dvi: main.tex \
	Figs/fig1.eps Figs/fig2.eps Figs/fig3.eps Figs/fig4.eps Figs/fig5.eps Figs/fig6.eps  Figs/fig7.eps
	platex main.tex

Figs/fig1.eps: Figs/setup.svgz
	inkscape -z -f Figs/setup.svgz -E Figs/fig1.eps
Figs/fig2.eps: Figs/rhod.svgz
	inkscape -z -f Figs/rhod.svgz -E Figs/fig2.eps
Figs/fig3.eps: Figs/Nyq_w17.svgz
	inkscape -z -f Figs/Nyq_w17.svgz -E Figs/fig3.eps
Figs/fig4.eps: Figs/z_typical.svgz
	inkscape -z -f Figs/z_typical.svgz -E Figs/fig4.eps
Figs/fig5.eps: Figs/circuit.svgz
	inkscape -z -f Figs/circuit.svgz -E Figs/fig5.eps
Figs/fig6.eps: Figs/T_CPE.svgz
	inkscape -z -f Figs/T_CPE.svgz -E Figs/fig6.eps
Figs/fig7.eps: Figs/pexp.svgz
	inkscape -z -f Figs/pexp.svgz -E Figs/fig7.eps
