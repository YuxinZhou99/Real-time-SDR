//
// Created by kaka on 3/19/21.
//

#include "dy4.h"
#include "rds.h"
#include "filter.h"
#include "logfunc.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>

float rds_fb_extr = 54e3;
float rds_fe_extr = 60e3;
float rds_fb_rec = 113.5e3;
float rds_fe_rec = 114.5e3;

float rds_Fs_mode0 = 2.4e5;
float rds_Fs_mode1 = 2.5e5;
int rds_bp_taps = 101;

//in order to get 57 kHz as ncoOutput
float rds_ncoScale = 0.5;
float rds_phaseAdjust = 0;
float rds_normBandwidth = 0.005;

int H[26][10] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
{0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
{1, 0, 1, 1, 0, 1, 1, 1, 0, 0},
{0, 1, 0, 1, 1, 0, 1, 1, 1, 0},
{0, 0, 1, 0, 1, 1, 0, 1, 1, 1},
{1, 0, 1, 0, 0, 0, 0, 1, 1, 1},
{1, 1, 1, 0, 0, 1, 1, 1, 1, 1},
{1, 1, 0, 0, 0, 1, 0, 0, 1, 1},
{1, 1, 0, 1, 0, 1, 0, 1, 0, 1},
{1, 1, 0, 1, 1, 1, 0, 1, 1, 0},
{0, 1, 1, 0, 1, 1, 1, 0, 1, 1},
{1, 0, 0, 0, 0, 0, 0, 0, 0, 1},
{1, 1, 1, 1, 0, 1, 1, 1, 0, 0},
{0, 1, 1, 1, 1, 0, 1, 1, 1, 0},
{0, 0, 1, 1, 1, 1, 0, 1, 1, 1},
{1, 0, 1, 0, 1, 0, 0, 1, 1, 1},
{1, 1, 1, 0, 0, 0, 1, 1, 1, 1},
{1, 1, 0, 0, 0, 1, 1, 0, 1, 1}};


void rds_fmPLL(std::vector<float> &ncoOut_I, std::vector<float> &ncoOut_Q, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst, float &triArg, float &trigOffset, float &ncoOut_buf_I, float &ncoOut_buf_Q){

    // scale factors for proportional/integrator terms
    float Cp = 2.666;
    float Ci = 3.555;

    //gain for the proportional term
    float Kp = normBandwidth * Cp;
    //gain for the integrator term
    float Ki = normBandwidth * normBandwidth * Ci;

    //initialize internal state
    ncoOut_I[0] = ncoOut_buf_I;
    ncoOut_Q[0] = ncoOut_buf_Q;

    for (auto i=0;i<pllIn.size();i+=1) {
        //phase detector
        float errorI = pllIn[i] * feedbackI; // complex conjugate of the
        float errorQ = pllIn[i] * (-1) * feedbackQ; //feedback complex exponential

        //four-quadrant arctangent discriminator for phase error detection
        float errorD = atan2(errorQ,errorI);

        //loop filter
        integrator += Ki * errorD;
        //update phase estimate
        phaseEst += (Kp * errorD) + integrator;
        //internal oscillator
        triArg = 2 * PI * (freq/Fs) * (trigOffset+i+1) + phaseEst;
        feedbackI = cos(triArg);
        feedbackQ = sin(triArg);
        ncoOut_I[i+1] = cos(triArg * ncoScale + phaseAdjust);
        ncoOut_Q[i+1] = sin(triArg * ncoScale + phaseAdjust);
    }
    trigOffset += pllIn.size();
    ncoOut_buf_I = ncoOut_I[ncoOut_I.size()-1];
    ncoOut_buf_Q = ncoOut_Q[ncoOut_Q.size()-1];

}

void bandpassRDS(std::vector<float> &h_extr, int mode){

    float Fs = 2.4e5;

    float normCenter = (rds_fe_extr+rds_fb_extr)/Fs;
    float normPass = (rds_fe_extr-rds_fb_extr)*2/Fs;

    for (auto i=0;i<rds_bp_taps;i++){
        if (i == (rds_bp_taps-1)/2){
            h_extr[i] = normPass;
        }
        else{
            h_extr[i] = normPass * (sin(PI*(normPass/2)*(i-(rds_bp_taps-1)/2)))/(PI*(normPass/2)*(i-(rds_bp_taps-1)/2));
        }

        h_extr[i] *= cos(PI*i*normCenter);
        h_extr[i] *= sin(i*PI/rds_bp_taps)*sin(i*PI/rds_bp_taps);
    }
}

void bpRDSRec(std::vector<float> &h_rec, int mode){

    float Fs = 2.4e5;

    float normCenter = (rds_fe_rec+rds_fb_rec)/Fs;
    float normPass = (rds_fe_rec-rds_fb_rec)*2/Fs;

    for (auto i = 0;i < rds_bp_taps;i += 1){
        if (i == (rds_bp_taps-1)/2){
            h_rec[i] = normPass;
        }
        else{
            h_rec[i] = normPass * (sin(PI*(normPass/2)*(i-(rds_bp_taps-1)/2)))/(PI*(normPass/2)*(i-(rds_bp_taps-1)/2));
        }

        h_rec[i] *= cos(PI*i*normCenter);
        h_rec[i] *= sin(i*PI/rds_bp_taps)*sin(i*PI/rds_bp_taps);
    }
}

void rds_mixer(std::vector<float> &rds_samples_I, std::vector<float> &rds_samples_Q, std::vector<float> &filter_demod_extr, std::vector<float> &ncoOut_I, std::vector<float> &ncoOut_Q){
    for (auto i=0; i<filter_demod_extr.size();i+=1){
        rds_samples_I[i] = filter_demod_extr[i] * ncoOut_I[i];
        rds_samples_Q[i] = filter_demod_extr[i] * ncoOut_Q[i];
    }
}

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1){
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

void rdsDataRecovery(std::vector<float> &rrc_outputs_I_3blocks, std::vector<float> &rrc_outputs_Q_3blocks, std::vector<float> &rrc_rec_I, std::vector<float> &rrc_rec_Q, int block_id, int &max_index){

    for (auto j=0;j<rrc_rec_I.size();j+=1) {
        rrc_rec_I[j] = rrc_outputs_I_3blocks[24 * j + max_index];
        rrc_rec_Q[j] = rrc_outputs_Q_3blocks[24 * j + max_index];
    }
}

void manchester(std::vector<int> &manchester_decoded_I, std::vector<float> &rrc_rec_I, int &buf_mach, int block_id){

    std::vector<int> manchester_I(38, 0);

    int count = 0;
    for (auto i=0;i<rrc_rec_I.size();i+=2){
        if (rrc_rec_I[i] < rrc_rec_I[i+1]) manchester_I[count] = 0;
        else if (rrc_rec_I[i] > rrc_rec_I[i+1]) manchester_I[count] = 1;
        count++;
    }

    //using saving state
    if (manchester_I[0] == buf_mach) manchester_decoded_I[0] = 0;
    else manchester_decoded_I[0] = 1;

    for (auto i=1;i<manchester_I.size();i+=1){

        if (manchester_I[i] == manchester_I[i-1]) manchester_decoded_I[i] = 0;
        else manchester_decoded_I[i] = 1;

        if (i == manchester_I.size()-1) buf_mach = manchester_I[i];
    }
}

void frame_sync(std::vector<int> &sync_I, int block_id, int &type_base, int &type_count, int &buf_type){

    int ini = 0;

    if (block_id == 3) ini = 26;
    else ini = 38;

    //  3     6      9     12       15     18
    //26-37 38-77 78-113 114-151 152-189 190-227
    for(auto i=ini;i<76;i++){
        std::vector<int> shift_I(26, 0);
        std::vector<int> offset(10, 0);
        int count = 0;
        //pick 26 samples from a whole bitstream
        //0 12 24...
        for(auto j=i-26;j<i;j++){
            shift_I[count] = sync_I[j];
            count++;
        }

        //matrix multiplication
        for (auto m= 0;m<10;m++){
            int single_result = 0;
            for (auto n=0;n<26;n++){
                single_result += shift_I[n] * H[n][m];
            }
            //determine if it is 0 or 1 for the syndrome results
            offset[m] = int(single_result%2);
        }

        if (offset[0] == 1 && offset[1] == 1 && offset[2] == 1 && offset[3] == 1 && offset[4] == 0 && offset[5] == 1 && offset[6] == 1 && offset[7] == 0 && offset[8] == 0 && offset[9] == 0) {
            if (type_base == 0 && type_count == 0) type_base = i%26;
            if (block_id == 3 || block_id == 6){
                if ((i-type_base)%26 == 0) {
                    std::cerr << "TYPE A " << i << std::endl;
                    buf_type = 0;
                }
                else std::cerr << "Probable false TYPE A " << i << std::endl;
            }
            else {
                if ((i+38*(block_id/3-2)-type_base)%26 == 0) {
                    std::cerr << "TYPE A " << i+38*(block_id/3-2) << std::endl;
                    buf_type = 0;
                }
                else std::cerr << "Probable false TYPE A " << i+38*(block_id/3-2) << std::endl;
            }
            type_count++;
        }
        if (offset[0] == 1 && offset[1] == 1 && offset[2] == 1 && offset[3] == 1 && offset[4] == 0 && offset[5] == 1 && offset[6] == 0 && offset[7] == 1 && offset[8] == 0 && offset[9] == 0) {
            if (type_base == 0 && type_count == 0) type_base = i%26;
            if (block_id == 3 || block_id == 6){
                if ((i-type_base)%26 == 0) {
                    std::cerr << "TYPE B " << i << std::endl;
                    buf_type = 1;
                }
                else std::cerr << "Probable false TYPE B " << i << std::endl;
            }
            else {
                if ((i+38*(block_id/3-2)-type_base)%26 == 0) {
                    std::cerr << "TYPE B " << i+38*(block_id/3-2) << std::endl;
                    buf_type = 1;
                }
                else std::cerr << "Probable false TYPE B " << i+38*(block_id/3-2) << std::endl;
            }
            type_count++;
        }
        if (offset[0] == 1 && offset[1] == 0 && offset[2] == 0 && offset[3] == 1 && offset[4] == 0 && offset[5] == 1 && offset[6] == 1 && offset[7] == 1 && offset[8] == 0 && offset[9] == 0 && buf_type != 2) {
            if (type_base == 0 && type_count == 0) type_base = i%26;
            if (block_id == 3 || block_id == 6){
                if ((i-type_base)%26 == 0) {
                    std::cerr << "TYPE C " << i << std::endl;
                    buf_type = 2;
                }
                else std::cerr << "Probable false TYPE C " << i << std::endl;
            }
            else {
                if ((i+38*(block_id/3-2)-type_base)%26 == 0) {
                    std::cerr << "TYPE C " << i+38*(block_id/3-2) << std::endl;
                    buf_type = 2;
                }
                else std::cerr << "Probable false TYPE C " << i+38*(block_id/3-2) << std::endl;
            }
            type_count++;
        }
        if (offset[0] == 1 && offset[1] == 1 && offset[2] == 1 && offset[3] == 1 && offset[4] == 0 && offset[5] == 0 && offset[6] == 1 && offset[7] == 1 && offset[8] == 0 && offset[9] == 0 && buf_type != 2) {
            if (type_base == 0 && type_count == 0) type_base = i%26;
            if (block_id == 3 || block_id == 6){
                if ((i-type_base)%26 == 0) {
                    std::cerr << "TYPE C' " << i << std::endl;
                    buf_type = 2;
                }
                else std::cerr << "Probable false TYPE C' " << i << std::endl;
            }
            else {
                if ((i+38*(block_id/3-2)-type_base)%26 == 0) {
                    std::cerr << "TYPE C' " << i+38*(block_id/3-2) << std::endl;
                    buf_type = 2;
                }
                else std::cerr << "Probable false TYPE C' " << i+38*(block_id/3-2) << std::endl;
            }
            type_count++;
        }
        if (offset[0] == 1 && offset[1] == 0 && offset[2] == 0 && offset[3] == 1 && offset[4] == 0 && offset[5] == 1 && offset[6] == 1 && offset[7] == 0 && offset[8] == 0 && offset[9] == 0) {
            if (type_base == 0 && type_count == 0) type_base = i%26;
            if (block_id == 3 || block_id == 6){
                if ((i-type_base)%26 == 0) {
                    std::cerr << "TYPE D " << i << std::endl;
                    buf_type = 3;
                }
                else std::cerr << "Probable false TYPE D " << i << std::endl;
            }
            else {
                if ((i+38*(block_id/3-2)-type_base)%26 == 0) {
                    std::cerr << "TYPE D " << i+38*(block_id/3-2) << std::endl;
                    buf_type = 3;
                }
                else std::cerr << "Probable false TYPE D " << i+38*(block_id/3-2) << std::endl;
            }
            type_count++;
        }
    }
}

//put in IF data as fm_demod
void rds(std::vector<int> &sync_I, std::vector<int> &sync_Q, std::vector<float> &rrc_outputs_I_3blocks, std::vector<float> &rrc_outputs_Q_3blocks, std::vector<float> &rrc_coeff, std::vector<float> &extr_coeff, std::vector<float> &fm_demod_rdsRec_coeff, std::vector<float> &rds_3k, std::vector<float> &rds_rationalLPF_coeff, int block_id, std::vector<float> &rds_samples_I, std::vector<float> &rds_samples_Q, std::vector<float> &fm_demod, std::vector<float> &rds_buf, std::vector<float> &rds_carrier_buf_I, std::vector<float> &rds_carrier_buf_Q, std::vector<float> &rds_carrier_buf, std::vector<float> &rds_up19_buf_I, std::vector<float> &rds_up19_buf_Q, std::vector<float> &rrc_buf_I, std::vector<float> &rrc_buf_Q, int block_size, float freq, int mode, float &feedbackI, float &feedbackQ, float &integrator, float &phaseEst, float &triArg, float &trigOffset, float &ncoOut_buf_I, float &ncoOut_buf_Q, int &max_index, int &buf_mach, int &type_base, int &type_count, int &buf_type){

    std::vector<float> filter_demod_extr(fm_demod.size(), 0.0);
    std::vector<float> filter_demod_rec(fm_demod.size(), 0.0);
    std::vector<float> squaring_results(fm_demod.size(), 0.0);
    std::vector<float> ncoOut_I(fm_demod.size()+1, 0.0);
    std::vector<float> ncoOut_Q(fm_demod.size()+1, 0.0);
    float Fs = 2.4e5;

    //RDS Extraction Channel
    stereoBPConv(filter_demod_extr, rds_buf, fm_demod, extr_coeff, fm_demod.size());//same as stereo
    int buf_count = 0;
    for (auto i = (fm_demod.size() - extr_coeff.size() + 1); i < fm_demod.size(); i+=1) {
        rds_buf[buf_count] = fm_demod[i];
        buf_count += 1;
    }

    //RDS Carrier Recovery
    //squaring non-linear is combined into conv process
    rdsRecBPFConv(filter_demod_rec, rds_carrier_buf, filter_demod_extr, fm_demod_rdsRec_coeff, squaring_results.size());//same as stereo
    buf_count = 0;
    for (auto i = (filter_demod_extr.size() - fm_demod_rdsRec_coeff.size() + 1); i < filter_demod_extr.size(); i+=1) {
        rds_carrier_buf[buf_count] = filter_demod_extr[i];
        buf_count += 1;
    }

    //pll
    rds_fmPLL(ncoOut_I, ncoOut_Q, filter_demod_rec, freq, Fs, rds_ncoScale, rds_phaseAdjust, rds_normBandwidth, feedbackI, feedbackQ, integrator, phaseEst, triArg, trigOffset, ncoOut_buf_I, ncoOut_buf_Q);

    //mixer
    rds_mixer(rds_samples_I, rds_samples_Q, filter_demod_extr, ncoOut_I, ncoOut_Q);

    //3k LPF combined with upsampled by 19
    std::vector<float> rds_3k_results_I_up19(rds_samples_I.size()*19, 0.0);
    std::vector<float> rds_3k_results_Q_up19(rds_samples_Q.size()*19, 0.0);
    rdsLPFConv_3k(rds_3k_results_I_up19, rds_carrier_buf_I, rds_samples_I, rds_3k, rds_samples_I.size());
    rdsLPFConv_3k(rds_3k_results_Q_up19, rds_carrier_buf_Q, rds_samples_Q, rds_3k, rds_samples_I.size());

    buf_count = 0;
    for (auto i=(rds_samples_I.size()-rds_3k.size()+1);i<rds_samples_I.size();i++){
        rds_carrier_buf_I[buf_count] = rds_samples_I[i];
        rds_carrier_buf_Q[buf_count] = rds_samples_Q[i];
        buf_count++;
    }

    //16k LPF combined with downsampled by 80
    //non-zero values in filter_rds_up19 = 608
    std::vector<float> filter_rds_down80_I(rds_3k_results_I_up19.size()/80,0.0);
    std::vector<float> filter_rds_down80_Q(rds_3k_results_Q_up19.size()/80,0.0);

    rds_16k_conv(filter_rds_down80_I, rds_up19_buf_I, rds_3k_results_I_up19, rds_rationalLPF_coeff, rds_3k_results_I_up19.size());
    rds_16k_conv(filter_rds_down80_Q, rds_up19_buf_Q, rds_3k_results_Q_up19, rds_rationalLPF_coeff, rds_3k_results_Q_up19.size());
    buf_count = 0;
    for (auto i=(rds_3k_results_I_up19.size()-rds_rationalLPF_coeff.size() + 1);i<rds_3k_results_I_up19.size();i++){
        rds_up19_buf_I[buf_count] = rds_3k_results_I_up19[i];
        rds_up19_buf_Q[buf_count] = rds_3k_results_Q_up19[i];
        buf_count++;
    }

    //RRC
    std::vector<float> rrc_outputs_I(filter_rds_down80_I.size(), 0.0);
    std::vector<float> rrc_outputs_Q(filter_rds_down80_Q.size(), 0.0);
    stereoBPConv(rrc_outputs_I, rrc_buf_I, filter_rds_down80_I, rrc_coeff, filter_rds_down80_I.size());
    stereoBPConv(rrc_outputs_Q, rrc_buf_Q, filter_rds_down80_Q, rrc_coeff, filter_rds_down80_Q.size());
    buf_count = 0;
    for (auto i=(filter_rds_down80_I.size()-rrc_coeff.size()+1);i<filter_rds_down80_I.size();i++){
        rrc_buf_I[buf_count] = filter_rds_down80_I[i];
        rrc_buf_Q[buf_count] = filter_rds_down80_Q[i];
        buf_count++;
    }

    if (block_id != 0){
        for (auto i=0;i<rrc_outputs_I.size();i++){
            rrc_outputs_I_3blocks[(block_id-1)%3*rrc_outputs_I.size()+i] = rrc_outputs_I[i];
            rrc_outputs_Q_3blocks[(block_id-1)%3*rrc_outputs_Q.size()+i] = rrc_outputs_Q[i];
        }
        if (block_id == 1){
            int wrong_cases = 0;
            for (auto i=0;i<24;i++){
                int wrong_cases_compare = 0;
                for (auto j=1;j<24;j+=2){
                    if ((rrc_outputs_I[(j-1)*24+i] > 0 && rrc_outputs_I[j*24+i] > 0) || (rrc_outputs_I[(j-1)*24+i] < 0 && rrc_outputs_I[j*24+i] < 0)) wrong_cases_compare++;
                }
                if (i==0) wrong_cases = wrong_cases_compare;
                else{
                    if (wrong_cases_compare < wrong_cases) {
                        max_index = i;
                        wrong_cases = wrong_cases_compare;
                    }
                }
            }
        }
    }

    //recovery
    //store 3 blocks of data, and starts from block 0
    std::vector<float> rrc_rec_I(1824/24, 0.0);//608*3/24
    std::vector<float> rrc_rec_Q(1824/24, 0.0);
    //Manchester
    std::vector<int> manchester_decoded_I((1824/2)/24, 0);//608*3/24

    if (block_id%3 == 0 && block_id != 0){
        rdsDataRecovery(rrc_outputs_I_3blocks, rrc_outputs_Q_3blocks, rrc_rec_I, rrc_rec_Q, block_id, max_index);
        manchester(manchester_decoded_I, rrc_rec_I, buf_mach, block_id);
        //dynamically changing sync_I size and shift data for 38 indexes as (block_id-3)%6 == 0
        if (block_id == 3){
            for (auto i=0;i<manchester_decoded_I.size();i++){
                sync_I[i] = manchester_decoded_I[i];
            }
        }
        else if (block_id == 6){
            std::vector<int> temp(38, 0);
            for (auto i=0;i<temp.size();i++){
                temp[i] = sync_I[i];
            }
            sync_I.resize(76, 0);
            for (auto i=0;i<temp.size();i++){
                sync_I[i] = temp[i];
            }
            for (auto i=0;i<manchester_decoded_I.size();i++){
                sync_I[38+i] = manchester_decoded_I[i];
            }
        }
        else {
            for (auto i=0;i<38;i++){
                sync_I[i] = 0;
            }
            std::rotate(sync_I.begin(), sync_I.begin()+38, sync_I.end());
            for (auto i=0;i<manchester_decoded_I.size();i++){
                sync_I[38+i] = manchester_decoded_I[i];
            }
        }
        frame_sync(sync_I, block_id, type_base, type_count, buf_type);
    }

//    float df = 1;
//    std::vector<float> fre = arange<float> (0, rrc_outputs_I_3blocks.size(), df);
//
//    if (block_id == 15){
//        logVector("rrc_output", rrc_rec_I, rrc_rec_Q); // log only positive freq
//        logVector("rrc_output1", fre, rrc_outputs_I_3blocks);
//        std::cerr << "Run: gnuplot -e 'set terminal png size 1024,768' example.gnuplot > ../data/example.png\n";
//    }
}