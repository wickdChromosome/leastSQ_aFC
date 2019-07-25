import pandas as pd
import numpy as np
from lmfit import Parameters
from lmfit import minimize
import pdb
import logging
from sklearn.linear_model import LinearRegression


def linear_estimate(h1,h2,expr,cov_dataf):

        """
        Use a covariate matrix to make an initial guess at the aFCs
        """

        #take only the samples that are in the haplotype matrix - this is the order the 
        #columns are organized in
        expr=expr[h1.columns].values
        cov_dataf = cov_dataf[h1.columns]
        cov_dataf = cov_dataf.append(h1 + h2)
        added_cov = cov_dataf
        X = (added_cov) * 0.5

        #get the intercepts,etc out of the linear fit
        reg = LinearRegression().fit(X.transpose(), expr)
        sa = reg.coef_ #the last one is the coef for the variant
        C0 = reg.intercept_ #the error

            
            
        #now store the results in the same matrix, to be passed to the least square optimizer
        cov_dataf = cov_dataf.transpose().append(pd.DataFrame(sa, columns=['covar_coeff'], index=cov_dataf.index).transpose()).transpose()

        #cov_dataf = cov_dataf.transpose().append(pd.DataFrame(np.full((cov_dataf.shape[1],1), 0), columns=['covar_coeff'], index=cov_dataf.index).transpose()).transpose()


        cov_dataf = cov_dataf.append(pd.DataFrame(np.full((cov_dataf.shape[1],1), C0), columns=['C0'], index=cov_dataf.columns).transpose()) #append C0s

        #cov_dataf = cov_dataf.append(pd.DataFrame(np.full((cov_dataf.shape[1],1), 0), columns=['C0'], index=cov_dataf.columns).transpose()) #append C0s


        return cov_dataf

def nocovar_linear_estimate(h1,h2,expr):

        """
        Make an initial guess at the aFCs, using a linear fit
        """

        #take only the samples that are in the haplotype matrix - this is the order the 
        #columns are organized in
        expr=expr[h1.columns].values
        added_cov = h1 + h2
        X = (added_cov) * 0.5

        #do the linear regression
        reg = LinearRegression().fit(X.transpose(),expr)
        sa = reg.coef_ #each coeff is for a different variant
        C0 = reg.intercept_ #the error

        sad = {}

        varc = 0
        for var in h1.index:

            sad[var] = sa[varc]
            varc += 1
            
        #return the aFCs as well as the C0
        return [sad,C0]

    
def gene_expression_normalization(gene_expression_df):
        
        """
        If the option is passed, normalize gene expression matrix
        """
        #replace nans
        expr = gene_expression_df.replace('',0).astype(float)
        expr_working = expr.copy()
        #drop genes if they have a median value of less than 100 counts over all samples for calculating the scaling factor 
        expr_working.drop([col for col, val in expr_working.median().iteritems() if val < 100], axis=1, inplace=True)
        #reference column
        ref_col = np.median(expr_working.values, axis=0)
        #div each value by the median
        med_div = expr_working/ref_col
        #take median of that over samples
        sample_med = np.median(med_div, axis=1)
        #now scale everything
        result = expr.transpose()/sample_med

        return result.transpose()

    
def fcn2fit(params,h1,h2,variants,expr,covarlist ):
    
        """
        Do least squares optimization with a covariate matrix
        """
    
        #unpack covariate list
        [id_cov,cov_dataf]= covarlist

        log2_c0=params['C0'].value
        sa=[]
        wa = []

        for arrayvalue in variants:

            sa.append(params[str(arrayvalue)].value)
           
        #get covariates only from matrix
        k= cov_dataf.transpose()[id_cov].transpose()

        #now make an array of params for each id
        for arrayvalue in id_cov:

            wa.append(params[str(arrayvalue)].value)

        #get expressions for only the common samples
        expr=expr[h1.columns].values
        h1 = h1.values
        h2 = h2.values

        #calculate model
        dotprod1 = np.dot(np.array(sa) , np.array(h1))
        dotprod2 = np.dot(np.array(sa) , np.array(h2))

        model = np.log2(2**(dotprod1) + 2**(dotprod2)) + log2_c0 + np.dot(k[k.columns[:-1]].transpose().values , wa) #also add the covar matrix

        #expressions were log transformed outside function
        error = (model - expr)

        #return squared error
        return error**2
    
        

    
def nocovar_fcn2fit(params,h1,h2,variants,expr):
    
        """
        Do least squares optimization 
        """
        
        log2_c0=params['C0'].value
        sa=[]

        for arrayvalue in variants:

            sa.append(params[str(arrayvalue)].value)

        #get expressions for only the common samples
        expr=expr[h1.columns].values
        h1 = h1.values
        h2 = h2.values

        #calculate model
        dotprod1 = np.dot(np.array(sa) , np.array(h1))
        dotprod2 = np.dot(np.array(sa) , np.array(h2))

        model = np.log2(2**(dotprod1) + 2**(dotprod2)) + log2_c0 
        #take the log2 of the expr (np.log2(expr)) if the input expressions were not log transformed!!!
        error = model - expr 


        return error**2


                        
def nonlin_solve(haplotype0_df, haplotype1_df, eqtl_dataf, expressions_df ,useful_genes, cov_dataf, neednorm, outname, is_cov, needlog):

        """
        Initializes the afcs with a linear fit, then uses those as a basis for the nonlinear fit.
        """
             
        #add a new column to eqtl_dataf to store afc values       
        eqtl_dataf['log2_aFC'] = np.full(len(eqtl_dataf.index.values),np.nan) 
        eqtl_dataf['log2_aFC_error'] = np.full(len(eqtl_dataf.index.values),np.nan) 
        eqtl_dataf['linear_log2_aFC'] = np.full(len(eqtl_dataf.index.values),np.nan) 
        eqtl_dataf['linear_c0'] = np.full(len(eqtl_dataf.index.values),np.nan) 
        eqtl_dataf['log2_aFC_c0'] = np.full(len(eqtl_dataf.index.values),np.nan) 
        
        #only use columns where genes are useful
        expressions_df = expressions_df[useful_genes]
        #replace nans with zeros if they exist in dataset
        expressions_df = expressions_df.fillna(0).astype(float)
        
        
        #if the user hasn't already normalized the data with deseq2 or another method
        if int(neednorm) == 0:

            print("normalizing dataset")
            expressions_df = gene_expression_normalization(expressions_df)
        
        #if the user hasn't already log transformed the data
        if int(needlog) == 0:
            
            print("Log transforming dataset")
            expressions_df = np.log2(expressions_df + 1)
        

                
        genecount = 0
        sample_out = False
        
        #if there was a covariant dataframe as input
        if is_cov != 0:

            #make the index the ID column, but keep the samples
            cov_dataf.index = cov_dataf[cov_dataf.columns[0]]
            cov_dataf = cov_dataf[cov_dataf.columns[1:]]

            #get number of covariates in matrix as index 
            id_cov = cov_dataf.index

            #iterate over gene in dataset, without the genes to skip
            for gene in useful_genes:

                #get all the variants for that one gene
                variants = eqtl_dataf[eqtl_dataf['gene_id'] == gene].variant_id.values

                #get the haplotype entries for these variant ids
                haplo0_variants = haplotype0_df[haplotype0_df.index.isin(variants)]
                haplo1_variants = haplotype1_df[haplotype1_df.index.isin(variants)]

                #get the expressions for this gene
                thisgene_expressions = expressions_df[gene]

                #if there are nan values for an individual for any of these variants in the haplo matrix, drop those individuals
                noNAN_cols = list(set(pd.concat([haplo0_variants, haplo1_variants], axis=1).dropna(axis=1).columns.values))
                haplo0_variants = haplo0_variants[noNAN_cols]
                haplo1_variants = haplo1_variants[noNAN_cols]
                
                thisgene_expressions = thisgene_expressions[thisgene_expressions.index.isin(noNAN_cols)]
                    
                #check again if there are variants for this gene to be sure
                if len(haplo0_variants.index.values) == 0:

                    print("no variants for this gene: " + str(gene))
                    continue
                    
                #for now, use 0 for initialization
                #get the initialization values
                covar = linear_estimate(haplo0_variants.astype(float), haplo1_variants.astype(float), thisgene_expressions.astype(float), cov_dataf)
                
                
                #initalize parameters object
                params = Parameters()
                     
                for variant in haplo0_variants.index.values:
                    
                    eqtl_dataf.linear_log2_aFC[(variant == eqtl_dataf.variant_id) & (eqtl_dataf.gene_id == gene)] = covar['covar_coeff'][variant] #coefficient for that variant
                    
                eqtl_dataf.linear_c0[eqtl_dataf.gene_id == gene] = covar['covar_coeff']['C0'] #save c0
                
                #add coeffs
                for parameter in covar['covar_coeff'].index:
                    params.add( str(parameter), value=covar['covar_coeff'][parameter])
                #do the actual minimization 
                result = minimize(fcn2fit, params, args=(haplo0_variants.astype(float), haplo1_variants.astype(float), haplo0_variants.index.values, thisgene_expressions.astype(float), [id_cov,covar]),method='leastsq', xtol=0.01, ftol=0.01)
                    
                #save parameters and results in results df
                for parameter in result.params:

                    #if parameter is an s value, variant_order excludes variants with nan values
                    if parameter in haplo0_variants.index.values:

                        eqtl_dataf.log2_aFC[(parameter == eqtl_dataf.variant_id) & (eqtl_dataf.gene_id == gene)] = result.params[parameter].value
                        eqtl_dataf.log2_aFC_error[(parameter == eqtl_dataf.variant_id) & (eqtl_dataf.gene_id == gene)] = result.params[parameter].stderr

                    elif parameter == 'C0': #its C0
                        eqtl_dataf.log2_aFC_c0[eqtl_dataf.gene_id == gene] = result.params[parameter]

                #keep track of genes for notifying the user
                genecount += 1   

                print("Done with genes: " + str(genecount), end='\r')
                if genecount > 500 and sample_out == False:

                    print("New aFCs out")
                    sample_out = True
                    eqtl_dataf.to_csv(outname + "_prelim", sep='\t')

                    
        #there has to be a more elegant way to do this
        else:
            
            #iterate over gene in dataset, without the genes to skip
            for gene in useful_genes:

                #get all the variants for that one gene
                variants = eqtl_dataf[eqtl_dataf['gene_id'] == gene].variant_id.values

                #get the haplotype entries for these variant ids
                haplo0_variants = haplotype0_df[haplotype0_df.index.isin(variants)]
                haplo1_variants = haplotype1_df[haplotype1_df.index.isin(variants)]

                #get the expressions for this gene
                thisgene_expressions = expressions_df[gene]

                #if there are nan values for an individual for any of these variants in the haplo matrix, drop those individuals
                noNAN_cols = list(set(pd.concat([haplo0_variants, haplo1_variants], axis=1).dropna(axis=1).columns.values))
                haplo0_variants = haplo0_variants[noNAN_cols]
                haplo1_variants = haplo1_variants[noNAN_cols]
                thisgene_expressions = thisgene_expressions[thisgene_expressions.index.isin(noNAN_cols)]

                #check again if there are variants for this gene to be sure
                if len(haplo0_variants.index.values) == 0:

                    print("no variants for this gene: " + str(gene))
                    continue

                #get the initialization values
                [sa, c0] = nocovar_linear_estimate(haplo0_variants.astype(float), haplo1_variants.astype(float), thisgene_expressions.astype(float))
                eqtl_dataf.linear_c0[eqtl_dataf.gene_id == gene] = c0 #save c0

                #initalize parameters object
                params = Parameters()

                #load initialization values into the parameters object, as well as store them in the results df
                for variant in list(haplo0_variants.index):

                    try:
                        
                        eqtl_dataf.linear_log2_aFC[(variant == eqtl_dataf.variant_id) & (eqtl_dataf.gene_id == gene)] = sa[variant] #coefficient for that variant

                    except:
                        pdb.set_trace()
                        
                    #add coefficient as new initial value
                    params.add( str(variant), value=0)#sa[variant], min=-6.64, max=6.64)

                params.add('C0', value=0)#c0 )

                #do the actual minimization 
                result = minimize(nocovar_fcn2fit, params, args=(haplo0_variants.astype(float), haplo1_variants.astype(float), haplo0_variants.index.values, thisgene_expressions.astype(float)), method='leastsq', xtol=0.01, ftol=0.01)

                #save parameters and results in results df
                for parameter in result.params:

                    #if parameter is an s value, variant_order excludes variants with nan values
                    if parameter in haplo0_variants.index.values:

                        eqtl_dataf.log2_aFC[(parameter == eqtl_dataf.variant_id) & (eqtl_dataf.gene_id == gene)] = result.params[parameter].value
                        eqtl_dataf.log2_aFC_error[(parameter == eqtl_dataf.variant_id) & (eqtl_dataf.gene_id == gene)] = result.params[parameter].stderr

                    else: #its C0
                        eqtl_dataf.log2_aFC_c0[eqtl_dataf.gene_id == gene] = result.params[parameter]

                #keep track of genes for notifying the user
                genecount += 1   

                print("Done with genes: " + str(genecount), end='\r')
                if genecount > 500 and sample_out == False:

                    print("New aFCs out")
                    sample_out = True
                    eqtl_dataf.to_csv(outname + "_prelim", sep='\t')

        print("Genes processed: " + str(genecount)) 
        eqtl_dataf.to_csv(outname, sep='\t', index=False)
