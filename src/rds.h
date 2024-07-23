//
// Created by kaka on 3/19/21.
//

#ifndef PROJECT_GROUP46_THURSDAY_RDS_H
#define PROJECT_GROUP46_THURSDAY_RDS_H
#include <iostream>
#include <vector>

void rds_fmPLL(std::vector<float> &, std::vector<float> &, std::vector<float> &, float, float, float, float, float, float &, float &, float &, float &, float &, float &, float &, float &);
void bandpassRDS(std::vector<float> &, int);
void bpRDSRec(std::vector<float> &, int);
void rds_mixer(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &);
void rdsDataRecovery(std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, int, int &);
void manchester(std::vector<float> &, std::vector<float> &, int &, int);
void frame_sync(std::vector<int> &, int, int &, int &, int &);
void rds(std::vector<int> &, std::vector<int> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, int, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, std::vector<float> &, int, float, int, float &, float &, float &, float &, float &, float &, float &, float &, int &, int &, int &, int &, int &);


#endif //PROJECT_GROUP46_THURSDAY_RDS_H
