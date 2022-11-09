@file:DependsOn("net.maizegenetics:tassel:5.2.84")

import net.maizegenetics.phenotype.Phenotype
import net.maizegenetics.phenotype.PhenotypeBuilder

val phenotypes: List<Phenotype> = emptyList() // make list of phenotypes
val intersectedPhenotype = PhenotypeBuilder().fromPhenotypeList(phenotypes).intersectJoin().build()[0]
