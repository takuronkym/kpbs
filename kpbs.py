#!/usr/bin/env python3

import sys
import numpy as np
import allel
import argparse
import pandas as pd


class PopulationSet():
    '''
    A class associated with generating a list represent correspondence between samples and populations.
    '''

    # A function to make a list object of sample names for a population based on a list, which is written in a file at gevin path.
    @staticmethod
    def _read_ids(fpath):
        li =[]
        with open(fpath) as fin:
            for line in fin:
                line = line.rstrip()
                if line != "":
                    li.append(line)
        return li

    # This generates a list object containing sample lists for each of three populations.
    @staticmethod
    def _make_subpop_li(p1, p2, p3):
        p_li = [[],[],[]]
        for i, path in enumerate([p1, p2, p3]):
            subp_li = PopulationSet._read_ids(path)
            p_li[i]= subp_li
        return p_li

    # This returns a list object containing index of samples in the input vcf file. Indices of samples for each population should be binned into a sub-list.
    @staticmethod
    def get_subpop_index(samples, p1, p2, p3):
        p_ind_li = [[],[],[]]
        subpops_li = PopulationSet._make_subpop_li(p1, p2, p3)
        for i,li in enumerate(subpops_li):
            for item in li:
                if item in samples:
                    smpl_ind = np.where(samples==item)[0][0]
                    p_ind_li[i].append(smpl_ind)
                else:
                    raise Exception("sample name [{}] was not found in the vcf.".format(item))
        return p_ind_li


class PbsDataSet(dict):
    '''
    PbsDataSet object is designed to include numpy ndarray of following items:
        Chromosome/Contig/Scaffold ID, SNP position, SNP ID, SNP quality, Reference for SNP sites, Fst values between three populations as well as PBS values for the populations.
    A certain index for these arrays corresponding to information for a single SNP site. Therefore these ndarrays must have equal length.
    '''
    def __init__(self,
                 cntg_arr,
                 snp_pos,
                 snp_id,
                 snp_qual,
                 snp_ref,
                 num12,
                 num23,
                 num31,
                 den12,
                 den23,
                 den31,
                 ori_ind=None,
                 pbs1=None,
                 pbs2=None,
                 pbs3=None,
                 pval1=None,
                 pval2=None,
                 pval3=None):

        super().__init__()
        self.key_li = ["cntg_arr", 
                    "snp_pos", 
                    "snp_id", 
                    "snp_qual", 
                    "snp_ref", 
                    "num12", 
                    "num23", 
                    "num31", 
                    "den12", 
                    "den23", 
                    "den31",
                    "ori_ind"]
        if ori_ind is None:
            ds_len = len(cntg_arr)
            ori_ind = np.array([i for i in range(ds_len)])
        self.update({
            "cntg_arr" : cntg_arr,
            "snp_pos" : snp_pos,
            "snp_id" : snp_id,
            "snp_qual" : snp_qual,
            "snp_ref" : snp_ref,
            "ori_ind" : ori_ind,
            "num12" : num12,
            "num23" : num23,
            "num31" : num31,
            "den12" : den12,
            "den23" : den23,
            "den31" : den31})

    # A function to calculate Fst between two populations.
    @staticmethod
    def _get_fst_comp(g, subpops, method):
        if method=="wc":
            a, b, c = allel.weir_cockerham_fst(g, subpops)
            #fsts = (np.sum(a, axis=1) / (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))
            num, den = np.sum(a, axis=1) , (np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1))
        elif method=="h":
            ac1 = g.count_alleles(subpop=subpops[0])
            ac2 = g.count_alleles(subpop=subpops[1])
            num, den = allel.hudson_fst(ac1, ac2)
            #fsts = num / den
        #return np.clip(fsts, None, 0.99999)
        return num,den
    
    @staticmethod
    def _get_variable_site_indx(cntg_gt_arr, subpops_ind, comp_pairs):
        inv_indx = []
        allel_var_li = [[],[],[]]
        for i_ in range(3):
            indxs = subpops_ind[i_]
            subp_gt_arr = cntg_gt_arr[:,indxs,:]
            # each item of allel_var_li consists of a tuple of ndarray: (allels found in the site, counts for each respective alleles)
            allel_var_li[i_] = [ np.unique(x, return_counts=True) for x in subp_gt_arr]
        for num, (x,y,z) in enumerate(zip(*allel_var_li)):
            li_ = [list(a[0]) for a in [x,y,z]]
            for ind_li in comp_pairs:
                allel_var_1 = li_[ind_li[0]]
                allel_var_2 = li_[ind_li[1]]
                # step for removal of -1 genotypes
                allel_var_1 = [x for x in allel_var_1 if x != (-1)]
                allel_var_2 = [x for x in allel_var_2 if x != (-1)]
                if (allel_var_1 == []) or (allel_var_2 == []): 
                    inv_indx.append(num)
                    break
                elif allel_var_1 == allel_var_2:
                    if len(allel_var_1) == 1:
                        inv_indx.append(num)
                        break
        var_indx = [ x for x in range(len(cntg_gt_arr)) if not x in inv_indx]
        return var_indx

    @staticmethod
    def _shuffled(arr):
        f_arr = arr.flatten()
        np.random.shuffle(f_arr)
        return f_arr.reshape([-1,2])

    @staticmethod
    def generate_dataset(callset, contig_id, subpops_ind, type_fst, shfl=None):
        '''
        This function generates a dataset object.
        This requires:
            - Callset data, which can be obtained by `allel.read_vcf()`
            - Contig/Chromosome/Scaffold ID to be input
            - a list including index information of three populations obtained with `PopulationSet.get_subpop_index()`.
        '''
        gt_data = callset["calldata/GT"]
        cntg_bools = (callset["variants/CHROM"] == contig_id)
        cntg_gt_data = gt_data[cntg_bools]
        if shfl=="shfl":
            cntg_gt_data = np.array([PbsDataSet._shuffled(x) for x in cntg_gt_data])
        cntg_gt_arr = allel.GenotypeArray(cntg_gt_data)
        
        # This list describes pairs of populations whose Fst is to be calculated.
        comp_pairs = [[0,1],
                      [1,2],
                      [2,0]]
        var_indx = PbsDataSet._get_variable_site_indx(cntg_gt_arr, subpops_ind, comp_pairs)
        num_d = {}
        den_d = {}
        for pair, k in zip(comp_pairs,["fst12","fst23","fst31"]):
            subpops = [subpops_ind[pair[0]], subpops_ind[pair[1]]]
            num, den = PbsDataSet._get_fst_comp(cntg_gt_arr, subpops, type_fst)
            num_d[k] = num
            den_d[k] = den
        ds = PbsDataSet(callset["variants/CHROM"][cntg_bools],
                          callset["variants/POS"][cntg_bools],
                          callset["variants/ID"][cntg_bools],
                          callset["variants/QUAL"][cntg_bools],
                          callset["variants/REF"][cntg_bools],
                          num_d["fst12"],
                          num_d["fst23"],
                          num_d["fst31"],
                          den_d["fst12"],
                          den_d["fst23"],
                          den_d["fst31"])
        # removal of sites with invariable allels pairs.
        for k in ds.key_li:
            ds[k] = ds[k][var_indx]
        return ds

    # This checks if any of three Fst for a SNP position is np.nan and then returns a dataset without such SNP position(s).
    def check_zero_fst_comp(self):
        for label, comps_arr in zip(["num","den"],
                        [[ self["num12"], self["num23"], self["num31"] ],
                         [ self["den12"], self["den23"], self["den31"] ]]):
            fst_arr = np.array(comps_arr)
            zero_ind = np.unique( np.where(fst_arr == 0)[1] )
            print("The number of sites with {0}=0: {1}".format(label, zero_ind.shape[0]),file=sys.stderr)
        '''
        args = []
        for k in self.key_li:
            args.append(np.delete(self[k],non_zero_ind))
        return PbsDataSet(*args)
        '''
    @staticmethod
    def _calc_pbs(fst_d):
        comp_sets = [["fst12","fst31","fst23"],
                     ["fst23","fst12","fst31"],
                     ["fst31","fst23","fst12"]]
        pbs_d = {}
        for pbs_id, comp_set in zip(["pbs1","pbs2","pbs3"], comp_sets):
            pbs_d[pbs_id] =  ((- np.log(1 - fst_d[comp_set[0]]) - np.log(1 - fst_d[comp_set[1]]) + np.log(1- fst_d[comp_set[2]]))/2)
        return pbs_d
    
    @staticmethod
    def _get_ave(pbs_arr,method="mean"):
        if method == "mean":
            return np.mean(pbs_arr)

    # A method to conduct window analysis of PBS values
    def get_window_pbs(self,window_size,step_size,start_pos=1):
        snp_pos = self["snp_pos"]
        window_size = window_size
        step_size = step_size

        s_pos = start_pos
        window_id = 1
        # a list to store lists to be input to PbsWindow(). The argments should be formatted in the following order: cntg_id, stt_pos, end_pos, wndw_id, snp_num, pbs1, pbs2, pbs3. All arguments are lists.
        pbs_window_input = [[] for _ in range(8)]
        while True:
            if s_pos <= snp_pos[-1]:
                contig_id = self["cntg_arr"][0]
                e_pos = s_pos + window_size - 1
                window_bools = np.where((snp_pos >= s_pos) & (snp_pos <= e_pos))
                snp_num = len(snp_pos[window_bools])
                sub_li = [contig_id, s_pos, e_pos, window_id, snp_num]
                fst_d = {}
                for num_k,den_k,fst_k in zip(["num12", "num23", "num31"], 
                                       ["den12", "den23", "den31"], 
                                       ["fst12", "fst23", "fst31"]):
                    pop_num_sub = self[num_k][window_bools]
                    pop_den_sub = self[den_k][window_bools]
                    if snp_num:
                        num_ave = PbsDataSet._get_ave(pop_num_sub)
                        den_ave = PbsDataSet._get_ave(pop_den_sub)
                        fst_d[fst_k] = np.clip(num_ave/den_ave, None, 0.99999)
                    else:
                        sub_li.append('NA')
                if snp_num:
                    pbs_d = PbsDataSet._calc_pbs(fst_d)
                    sub_li += [pbs_d[k] for k in ["pbs1","pbs2","pbs3"]]
                    del pbs_d
                else:
                    pass
                for ind, item in enumerate(sub_li):
                    pbs_window_input[ind].append( sub_li[ind] )
                s_pos += step_size
                window_id += 1
            else:
                break
        wndw_ds = PbsWindow(*pbs_window_input)
        return wndw_ds
    
    def print_pbs_table(self):
        head_li = ["Chromosome/Contig ID", "position", "ID", "PBS for population1", "PBS for population2", "PBS for population3", "P-value for pop1 PBS", "P-value for pop2 PBS", "P-value for pop3 PBS"]
        print( "\t".join(head_li) )
        key_li = ["cntg_arr", "snp_pos", "snp_id", "pbs1", "pbs2", "pbs3", "pval1", "pval2", "pval3"]
        for i in range( len(self["cntg_arr"]) ):
            out_li = []
            for k in key_li:
                out_li.append( "NA" if self[k][i] is None else self[k][i] )
            out_li = list(map(str, out_li))
            print( "\t".join(out_li) )

            
class PbsWindow(dict):
    def __init__(self,
                cntg_id,
                stt_pos,
                end_pos,
                wndw_id,
                snp_num,
                pbs1,
                pbs2,
                pbs3,
                pval1=None,
                pval2=None,
                pval3=None):
        super().__init__()
        self.key_li = ["cntg_id", "stt_pos", "end_pos", "wndw_id", "snp_num", "pbs1", "pbs2", "pbs3", "pval1", "pval2", "pval3"]
        self.update({"cntg_id" : cntg_id,
                    "wndw_id" : wndw_id,
                    "stt_pos" : stt_pos,
                    "end_pos" : end_pos,
                    "snp_num" : snp_num,
                    "pbs1" : pbs1,
                    "pbs2" : pbs2,
                    "pbs3" : pbs3,
                    "pval1": pval1,
                    "pval2": pval2,
                    "pval3": pval3})
    
    def output_wndw_table(self):
        table_len = len(self["wndw_id"])
        if self["pval1"] is None:
            for pval in ["pval1", "pval2", "pval3"]:
                self[pval] = ["NA" for _ in range(table_len)]
        out = []
        for ind in range( table_len ):
            out_li = []
            for k in self.key_li:
                out_li.append( self[k][ind] )
            out.append(out_li)
        return out

    def generate_condenced(self):
        pbs_window_input = [[] for _ in range(8)]
        table_len = len(self["wndw_id"])
        for ind in range(table_len):
            sub_li = []
            if self["pbs1"][ind] != "NA":
                for k in self.key_li[0:8]:
                    sub_li.append(self[k][ind])
                for ind, item in enumerate(sub_li):
                    pbs_window_input[ind].append( sub_li[ind] )
        cond_wndw_ds = PbsWindow(*pbs_window_input)
        return cond_wndw_ds

    # A method to calculate p-values based on the Monte-Carlo simulation
    def run_pbs_mc(self, callset, contig_id, subpops_ind, type_fst, num_itr):
        ori_len = int(  self["wndw_id"].pop() )
        mc_data = PbsMC(contig_id, ori_len, num_itr)
        [mc_data._update_window(callset, contig_id, subpops_ind, type_fst, num) for num in range(num_itr)]
        mc_data._finalize()

        ori_ind_li = [(num - 1) for num in self["wndw_id"]]
        itr_li =[[self["pbs1"], mc_data.mcp1_li, "pval1"],
                 [self["pbs2"], mc_data.mcp2_li, "pval2"], 
                 [self["pbs3"], mc_data.mcp3_li, "pval3"]]
        for ind in range(3):
            [pbs_li, mc_li, p_key] = [*itr_li[ind]]
            tmp_pval_li = []
            for ori_ind, pbs in zip(ori_ind_li, pbs_li):
                bool_li = ~np.isnan(mc_li[ori_ind])
                nonan_mc_li = mc_li[ori_ind][bool_li]
                # counting numbers of the number of items that are greater or equal to the actual pbs value
                num_item =  np.where(nonan_mc_li >= pbs)[0].shape[0] 
                # calculating a percentile for the range from a maximum value in the simulated pbs to a actual pbs value. 1 is added to account the actual pbs
                pval = ( (num_item + 1 ) / (mc_data.sum_li[ind][ori_ind] + 1 ) )
                tmp_pval_li.append(pval)
            self[p_key] = np.array(tmp_pval_li)


class Output():
    def __init__(self):
        self.head_li = ["#Chromosome/Contig ID", "start position of window", "end position of window", "window ID", "number of SNP in the window", "PBS for population1", "PBS for population2", "PBS for population3", "P-value for pop1 PBS", "P-value for pop2 PBS", "P-value for pop3 PBS"]
        self.out_li = []
    def add_line(self, out_li):
        self.out_li += out_li
    def create_df(self):
        self.out_df = pd.DataFrame(data=self.out_li, columns=self.head_li)
    def add_bs(self, bs_pbs_d):
        for pbs_column, pvalue_column, bs_k in zip(["PBS for population1", "PBS for population2", "PBS for population3"],
                                                   ["P-value for pop1 PBS", "P-value for pop2 PBS", "P-value for pop3 PBS"],
                                                   ["pbs1", "pbs2", "pbs3"]):
            self.out_df[pvalue_column] = ['NA' if pbs == 'NA' 
                                          else np.count_nonzero(bs_pbs_d[bs_k] >= pbs)/len(bs_pbs_d[bs_k]) 
                                          for pbs in self.out_df[pbs_column].values]
    def output(self):
        self.out_df.to_csv(sys.stdout, index=False, sep="\t")


class PbsMC():
    def __init__(self,contig_id,ori_len,num_itr):
        self.contig_id = contig_id
        self.mcp1_li =[np.full(num_itr,np.nan) for i in range(ori_len)]
        self.mcp2_li =[np.full(num_itr,np.nan) for i in range(ori_len)]
        self.mcp3_li =[np.full(num_itr,np.nan) for i in range(ori_len)]
        self.sum_li = [[],[],[]]
    
    def _append_mcpbs(self,ori_ind, num, mcp1,mcp2,mcp3):
        for mcp_li, mcp in zip([self.mcp1_li, self.mcp2_li, self.mcp3_li],[mcp1, mcp2, mcp3]):
            mcp_li_len = len(mcp_li)
            for ind,pbs in zip(list(ori_ind),mcp):
                if ind < mcp_li_len:
                    mcp_li[ind][num] = pbs.astype("float32")

    def _update_window(self, callset, contig_id, subpops_ind, type_fst, num):
        ds = PbsDataSet.generate_dataset(callset, contig_id, subpops_ind, type_fst, "shfl")
        wndw_ds = ds.get_window_pbs(window_size, step_size)
        cond_wndw_ds = wndw_ds.generate_condenced()
        ori_ind = [(num - 1) for num in cond_wndw_ds["wndw_id"]]
        mcp1, mcp2, mcp3 = cond_wndw_ds["pbs1"], cond_wndw_ds["pbs2"], cond_wndw_ds["pbs3"]
        self._append_mcpbs(ori_ind, num, mcp1, mcp2, mcp3)
        
    def _finalize(self):
        for i,mcp in enumerate([self.mcp1_li, self.mcp2_li, self.mcp3_li]):
            for item in mcp:
                bool_li = ~np.isnan(item)
                nonan = item[bool_li]
                self.sum_li[i].append(nonan.shape[0])


class Bootstrap():
    def __init__(self, gff_path):
        self.gff_df = allel.gff3_to_dataframe(gff_path)
        self.numden = {k:[[], #num
                          []] #den
                       for k in ['num12den12', 'num23den23', 'num31den31']}
    def intergenic(self, contig, ds):
        df = self.gff_df[self.gff_df["type"] == "gene"]
        if len(df) == 0:
            raise Exception
        else:
            pass
        df = df[df["seqid"] == contig]
        if len(df) == 0:
            raise Exception
        else:
            pass
        gene_regions = df[["start", "end"]].values
        gene_regions[:,0] -= (5000 - 1)
        gene_regions[:,1] += (5000 - 1)
        
        snp_pos = ds['snp_pos']
        bools = [(snp_pos >= region[0])&(snp_pos <= region[1]) for region in gene_regions]
        self.bools = np.all(~np.array(bools), axis=0)
    def update_(self, ds):
        for num,den in zip(['num12', 'num23', 'num31'], ['den12', 'den23', 'den31']):
            num_li = ds[num][self.bools]
            den_li = ds[den][self.bools]
            self.numden[num+den] = np.hstack([self.numden[num+den], [num_li,den_li]])
    def bootstrap(self, num_bs, window_size):
        index = self.numden.values().__iter__().__next__().shape[1]
        self.bs_pbs_d = {k:np.array([]) for k in ['pbs1', 'pbs2', 'pbs3']}
        for _ in range(num_bs):
            bs_index = np.random.choice(range(index), window_size, replace=True)
            bs_fst_d = {}
            for k, fst_k in zip(['num12den12', 'num23den23', 'num31den31'], ['fst12', 'fst23', 'fst31']):
                bs_numden = self.numden[k][:, [bs_index]]
                bs_num_ave = PbsDataSet._get_ave(bs_numden[0])
                bs_den_ave = PbsDataSet._get_ave(bs_numden[1])
                bs_fst_d[fst_k] = np.clip(bs_num_ave/bs_den_ave, None, 0.99999)
            bs_pbs_d_sub = PbsDataSet._calc_pbs(bs_fst_d)
            for k in ['pbs1', 'pbs2', 'pbs3']:
                self.bs_pbs_d[k] = np.append(self.bs_pbs_d[k], bs_pbs_d_sub[k])

def main():
    parser = argparse.ArgumentParser(description='Script to calculate PBS between specified three populations from an input vcf file') 
    parser.add_argument('--vcf', help='vcf file')
    parser.add_argument('--pop1', help='sample ID list for population1')
    parser.add_argument('--pop2', help='sample ID list for population2')
    parser.add_argument('--pop3', help='sample ID list for population3')
    parser.add_argument('-w','--window_size',type=int, help='window size in nucleotide length (bp) to be used in the window analysis')
    parser.add_argument('-s','--step_size',type=int, help='step size in nucleotide length (bp) to be used in the window analysis')
    #parser.add_argument('-t','--type_fst',type=str, choices=['wc', 'h'], default='h', help='type of Fst estimator (optional); Weir & Cockerham ("wc") or Hudson ("h"; default)')
    parser.add_argument('-mc','--num_mc_simulation',type=int, help='the number of sampling in the allele-shuffiling Monte-Carlo simulation')
    parser.add_argument('--gff', help='GFF file')
    parser.add_argument('--num_bs', type=int, help='Number of bootstrap sampling')

    args = parser.parse_args()

    vcf_path = args.vcf
    subpop1_path = args.pop1
    subpop2_path = args.pop2
    subpop3_path = args.pop3
    global window_size
    window_size = args.window_size
    global step_size
    step_size = args.step_size
    #type_fst = args.type_fst
    type_fst = "h"
    num_mc = args.num_mc_simulation
    gff_path = args.gff
    num_bs = args.num_bs

    # a global variable for a future implementation of parallel processing
    global thread_num
    thread_num = 1


    if step_size > window_size:
        raise Exception("The step size must be smaller than or equal to the window size.\ncurrent setting: step_size={0} window_size={1}".format(step_size,window_size))
    if step_size <= 0 or window_size <= 0:
        raise Exception("The step size and the window size must be integers greater than 0.")
    
    if num_bs is not None:
        bootstrap = Bootstrap(gff_path)

    callset = allel.read_vcf(vcf_path,
                             fields=['samples',
                             'variants/CHROM', 
                             'variants/POS', 
                             'variants/ID', 
                             'variants/QUAL', 
                             'variants/REF', 
                             'calldata/GT'])
    contig_li = list(callset['variants/CHROM'])
    contig_arr = np.array(sorted(set(contig_li), key=contig_li.index))
    del contig_li
    output = Output()
    for contig in contig_arr:
        subpops_ind = PopulationSet.get_subpop_index(callset["samples"], subpop1_path, subpop2_path, subpop3_path)
        ds = PbsDataSet.generate_dataset(callset, contig, subpops_ind, type_fst)
        ds.check_zero_fst_comp()
        wndw_ds = ds.get_window_pbs(window_size, step_size)        
        if num_mc is not None:
            wndw_ds = wndw_ds.generate_condenced()
            wndw_ds.run_pbs_mc(callset, contig, subpops_ind, type_fst, num_mc)
        out_li = wndw_ds.output_wndw_table()
        output.add_line(out_li)
        if num_bs is not None:
            bootstrap.intergenic(contig, ds)
            bootstrap.update_(ds)
    output.create_df()
    if num_bs is not None:
        bootstrap.bootstrap(num_bs, window_size)
        output.add_bs(bootstrap.bs_pbs_d)

    output.output()

if __name__ == '__main__':
    main()
