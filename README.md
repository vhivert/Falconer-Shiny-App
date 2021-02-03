<h1>Welcome in the Falconer ShinyApp</h1><br><br>
The aim of this App is to show how the combination of gene action and allele frequencies at causal loci translate to genetic variance and genetic variance components for a complex trait. Although the theory underlying the App is more than a century old, it is highly relevant in the current era of genome-wide association studies (GWAS). The App can be used to demonstrate the relationship between a SNP effect size estimated from GWAS and the variation the SNP generates in the population, <i>i.e.</i> how locus-specific effects lead to individual differences. In addition, it can also be used to demonstrate how within and between locus interactions (dominance and epistasis, respectively) usually do not lead to a large amount of non-additive variance relative to additive variance, and therefore that these interactions usually do not explain individual differences in a population.<br><br>

The three models described below mainly illustrate the Chapters 7 and 8 of Falconer and Mackay (1996) and the Chapter 5 of Lynch and Walsh (1998).<br><br>

<h3><b>Single-locus Model with additive and dominance effect:</b></h3><br>

In this single-locus model, we consider a biallelic locus with allele A<sub>1</sub> and A<sub>2</sub> in frequencies <i>p</i> and 1-<i>p</i>. Under panmixia (<i>i.e.</i> random mating) and Hardy-Weinberg equilibrium, the expected genotype frequencies are (1-<i>p</i>)<sup>2</sup>,2<i>p</i> (1-<i>p</i>) and <i>p</i><sup>2</sup>, for A<sub>2</sub>A<sub>2</sub>, A<sub>1</sub>A<sub>2</sub> and A<sub>1</sub>A<sub>1</sub> respectively. 
      We arbitrarily assign genotypic values (the mean trait values for the different genotype class) -<i>a</i>, <i>d</i> and <i>a</i> to the three genotypes, <i>d</i> representing the dominance effect (within locus interaction, no interaction when <i>d</i> = 0) and 2<i>a</i> the difference between the two homozygotes. Under this model, the population mean is:<br><br>

M = (2<i>p</i>-1)<i>a</i> + 2<i>p</i>(1-<i>p</i>)<i>d</i><br><br>

<b>Average effect of gene (allele) substitution (also called additive effect in the literature)</b><br>

The transmission of value from parents to offspring occurs through their genes (alleles) and not their genotypes. The average effect of gene substitution (&#120572;) is defined as the average effect on the trait when substituting alleles at this locus in the population. It can also be defined as the mean value of genotypes produced by different gametes:<br><br>

&#120572; = <i>a</i> + (1-2<i>p</i>)<i>d</i> <br><br>

Importantly, &#120572; is also the slope of the linear regression of the genotype means, weighted by their frequency, on the A<sub>1</sub> allele dosage (0, 1 or 2). <br><br>

When performing a standard GWAS, individual phenotypes <i>y</i> are regressed on the number <i>x</i> (<i>x</i> = 0, 1, 2) of reference alleles at a given locus, , <i>i.e.</i>, the allelic “dosage”, where the reference allele for this dosage count is arbitrarily the major or minor allele (but this arbitrary choice is reflected in the sign of the regression coefficient &beta;:<br><br>

y = &mu; + &beta;x + e<br><br>

Where the residuals e include both the non-additive genetic effects at the locus, the genetic effects (additive and non-additive) at other loci and an environmental and/or chance effect (non-genetic). The quantity of interest is the slope &beta; of the model (the effect size of the locus), which is the average effect of allele substitution, hence &beta; = &#120572;.<br><br>

<b>Additive (breeding) values and dominance deviations</b><br>

The breeding values are the expected genotypic values under additivity (the predictions from the linear model). Expressed as deviations from the population mean M, the breeding values of the 3 genotypes A<sub>2</sub>A<sub>2</sub>, A<sub>1</sub>A<sub>2</sub> and A<sub>1</sub>A<sub>1</sub> are
      -2<i>p</i>&#120572;, (1-2<i>p</i>)&#120572; and (2-2<i>p</i>)&#120572;. The residuals of the linear regression are the deviations due to the within locus interaction (dominance)<br><br>

<b>Genetic variance</b><br>

The total genotypic variance (<i>V<sub>G</sub></i>) is partitioned into Additive (<i>V<sub>A</sub></i>) and Dominance (<i>V<sub>D</sub></i>) variance.<br><br>

<i>V<sub>G</sub></i> = <i>V<sub>A</sub></i> + <i>V<sub>D</sub></i><br><br>

<b>Additive variance</b><br>

The Additive variance (<i>V<sub>A</sub></i>) is the variance of additive (breeding) values. When values are expressed in terms of deviation from the population mean, the variance simply become the mean of the squared values. Hence, <i>V<sub>A</sub></i> is obtained by squaring the additive (breeding) values described above, multiplying by the corresponding frequencies and summing over the 3 genotypes, leading to:<br><br>

<i>V<sub>A</sub></i> = 2<i>p</i>(1-<i>p</i>)&#120572;<sup>2</sup> = 2<i>p</i>(1-<i>p</i>)[<i>a</i>+<i>d</i>(1-2<i>p</i>)]<sup>2</sup> = <i>H</i>&#120572;<sup>2</sup>, <br><br>

with <i>H</i> the heterozygosity at the locus. Note that <i>V<sub>A</sub></i> is the variance explained by a SNP in GWAS (2<i>p</i>(1-<i>p</i>)&#120572;<sup>2</sup> = 2<i>p</i>(1-<i>p</i>)&beta;<sup>2</sup>) and contain both a term due to additivity (<i>a</i>) and dominance (<i>d</i>) through the average effect &#120572;.<br><br>

<b>Dominance variance</b><br>

Similarly, the dominance variance (<i>V<sub>D</sub></i>) is the variance of dominance deviations:<br><br>

<i>V<sub>D</sub></i> = (2<i>p</i> (1-<i>p</i>)<i>d</i>)<sup>2</sup> = <i>H</i><sup>2</sup><i>d</i><sup>2</sup><br><br>

Therefore, the dominance variance disproportionally depends on the locus heterozygosity compared to the additive variance (<i>H</i><sup>2</sup> versus <i>H</i>).<br><br>

<h3><b>Two-locus Model with additive and additive-by-additive effect:</b></h3><br>

We extend the one-locus to a two-locus model with additive and additive-by-additive epistatic interaction only, assuming no within loci dominance effects (<i>d</i> = 0 at both loci). We introduce a second (unlinked) locus with alleles B<sub>1</sub> and B<sub>2</sub> in frequencies <i>q</i> and 1-<i>q</i> respectively. The genotypic values and allele frequencies of the 9 genotypes are:
      <br>
       <table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>A<sub>2</sub>A<sub>2</sub></th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>A<sub>1</sub>A<sub>2</sub></th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>A<sub>1</sub>A<sub>1</sub></th>
</tr>
</thead>
<tbody>
<tr>
<td style='border-right: 1px solid black; text-align: left;'><b>B<sub>2</sub>B<sub>2</sub></b></td>
<td style='text-align: center;'><b>-<i>a</i><sub>A</sub>-<i>a</i><sub>B</sub>+<i>a</i><sub>AB</sub></b><br>&nbsp;&nbsp;&nbsp;(1-<i>p</i>)<sup>2</sup>(1-<i>q</i>)<sup>2</sup>&nbsp;&nbsp;&nbsp;</td>
<td style='text-align: center;'><b>-<i>a</i><sub>B</sub></b><br>2<i>p</i>(1-<i>p</i>)(1-<i>q</i>)<sup>2</sup></td>
<td style='text-align: center;'><b><i>a</i><sub>A</sub>-<i>a</i><sub>B</sub>-<i>a</i><sub>AB</sub></b><br><i>p</i><sup>2</sup>(1-<i>q</i>)<sup>2</sup></td>
</tr>
<tr>
<td style='border-right: 1px solid black; text-align: left;'><b>B<sub>1</sub>B<sub>2</sub></b></td>
<td style='text-align: center;'><b>-<i>a</i><sub>A</sub></b><br>(1-<i>p</i>)<sup>2</sup>2<i>q</i>(1-<i>q</i>)</td>
<td style='text-align: center;'><b>0</b><br>4<i>p</i>(1-<i>p</i>)<i>q</i>(1-<i>q</i>)</td>
<td style='text-align: center;'><b><i>a</i><sub>A</sub></b><br>2<i>p</i><sup>2</sup><i>q</i>(1-<i>q</i>)</td>
</tr>
<tr>
<td style='border-bottom: 2px solid grey; border-right: 1px solid black; text-align: left;'><b>B<sub>1</sub>B<sub>1</sub></b></td>
<td style='border-bottom: 2px solid grey; text-align: center;'><b>-<i>a</i><sub>A</sub>+<i>a</i><sub>B</sub>-<i>a</i><sub>AB</sub></b><br>(1-<i>p</i>)<i>q</i><sup>2</sup></td>
<td style='border-bottom: 2px solid grey; text-align: center;'><b><i>a</i><sub>B</sub></b><br>2<i>p</i>(1-<i>p</i>)<i>q</i><sup>2</sup></td>
<td style='border-bottom: 2px solid grey; text-align: center;'><b><i>a</i><sub>A</sub>+<i>a</i><sub>B</sub>+<i>a</i><sub>AB</sub></b><br><i>p</i><sup>2</sup><i>q</i><sup>2</sup></td>
</tr>
</tbody>
</table>
where <i>a</i><sub>A</sub> (<i>a</i><sub>B</sub>) is the genotypic value for the upper homozygote A<sub>1</sub>A<sub>1</sub> (B<sub>1</sub>B<sub>1</sub>) and <i>a</i><sub>AB</sub> is the additive-by-additive interaction effect. This is a re-parametrization of the model described by Mäki-Tanila and Hill (2014).<br><br>

<b>Population mean</b><br>

In our model, the mean of the genotypic values is:<br><br>

M = <i>a</i><sub>A</sub>(2<i>p</i>-1) + <i>a</i><sub>B</sub>(2<i>q</i>-1) + <i>a</i><sub>AB</sub>(1-2(<i>p</i>+<i>q</i>)+4<i>pq</i>)<br><br>

Note that the expression of M depends on the arbitrarily assigned genotypic values.<br><br>

<b>Average effect of gene (allele) substitution</b><br>

In this model, the locus specific average effects are:<br><br>

&#120572;<sub>A</sub> = <i>a</i><sub>A</sub> + 2<i>qa</i><sub>AB</sub><br>

&#120572;<sub>B</sub> = <i>a</i><sub>B</sub> + 2<i>pa</i><sub>AB</sub><br><br>

<b>Genetic variance</b><br>

The total genotypic variance (<i>V<sub>G</sub></i>) of the model is partitioned into Additive (<i>V<sub>A</sub></i>) and Additive-by-Additive (<i>V<sub>AA</sub></i>) variance.<br><br>

<i>V<sub>G</sub></i>=<i>V<sub>A</sub></i> + <i>V<sub>AA</sub></i><br><br>

<b>Additive variance</b><br>

The additive variance of the model is:<br><br>

<i>V<sub>A</sub></i> = &sum;<sub>i</sub><i>H</i><sub>i</sub>&#120572;<sub>i</sub><sup>2</sup>, with <i>H</i><sub>i</sub> the heterozygosity at locus i and &#120572;<sub>i</sub> the average effect of locus i. Hence:<br><br>

<i>V<sub>A</sub></i> = 2<i>p</i>(1-<i>p</i>)[<i>a</i><sub>A</sub>+2<i>qa</i><sub>AB</sub>] + 2<i>q</i>(1-<i>q</i>)[<i>a</i><sub>B</sub>+2<i>pa</i><sub>AB</sub>] <br><br>

Note that <i>V<sub>A</sub></i> contains a term due pairwise additive-by-additive interaction between locus A and B (<i>a</i><sub>AB</sub>).<br><br>

<b>Additive-by-Additive variance</b><br>

The additive-by-additive variance of the model is:<br><br>

<i>V<sub>AA</sub></i> = &sum;<sub>i</sub>&sum;<sub>j>i</sub><i>H</i><sub>i</sub><i>H</i><sub>j</sub><i>a</i><sub>ij</sub><sup>2</sup>, with <i>H</i><sub>i</sub> the heterozygosity at locus i and <i>a</i><sub>ij</sub> the additive-by-additive interaction effect between locus i and j. Hence:<br><br>

<i>V<sub>AA</sub></i> = 4<i>p</i>(1-<i>p</i>)<i>q</i>(1-<i>q</i>)<i>a</i><sub>AB</sub><br><br>

Therefore, the additive-by-additive variance disproportionally depends on the locus heterozygosity as compared to the additive variance.<br><br>

<h3><b>General two locus model:</b></h3><br>

Lastly, we use a generalized two-locus model where the user can provide all the genotypic values in an interactive table and choose the allele frequencies at the two loci (<i>p</i> and <i>q</i>). The genotypic values as well as the linear regressions are plotted as a function of the A<sub>1</sub> allelic dosage for the different genotypes at locus B, as well as the linear regression of the genotypic values weighted by their frequency on the A<sub>1</sub> allele dosage. The total genotypic variance (<i>V<sub>G</sub></i>) of this model is then partitioned in five components:<br><br>

<i>V<sub>G</sub></i> = <i>V<sub>A</sub></i> + <i>V<sub>D</sub></i> + <i>V<sub>AA</sub></i> + <i>V<sub>AD</sub></i> + <i>V<sub>DD</sub></i><br><br>

Where <i>V<sub>AD</sub></i> is the additive-by-dominance variance and <i>V<sub>DD</sub></i> the dominance-by-dominance variance. We use the least square approach described in Lynch and Walsh (1998) Chapter 5 to derive the different variance components and display their values.
  