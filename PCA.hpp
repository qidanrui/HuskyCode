#include <algorithm>
#include <cfloat>
#include <string>
#include <vector>
#include <functional>

#include "core/context.hpp"
#include "core/engine.hpp"
#include "lib/dcaaggregator.hpp"
#include "lib/topk.hpp"
#include "gperftools/profiler.h"
#include "modules/inputformat/lineinputformat.hpp"

class PCAObject {
public:
  int label;
  std::vector<double> inputFeatures;
  
  typedef int KeyT;

  virtual KeyT const & id() const {
    return label;
  }
};

class covMatrixObject{
public:
  int label;
  std::vector<double> row;
  std::vector<double> col;
  std::vector<double> eigenvectorRow;
  std::vector<double> eigenvectorCol;

  double eigenValue;

  typedef int KeyT;
  virtual KeyT const & id() const {
    return label;
  }
};

class PCA {
private:
  int degree;
  int vectorNum;
  int inputVectorDegree;
  int nRow;
  int nCol;
  int nCount;
  int dbMax;

  double dbApp;  
  double dbApq;  
  double dbAqq;
  double dbAngle;  
  double dbSinTheta;  
  double dbCosTheta;  
  double dbSin2Theta;  
  double dbCos2Theta; 

  Husky::ObjList<PCAObject>& inputList;
  Husky::ObjList<covMatrixObject>& covMatrix;
  
public:
  PCA(Husky::ObjList<PCAObject>& inputList, Husky::ObjList<covMatrixObject>& covMatrix) : inputList(inputList), covMatrix(covMatrix){
    nRow = 0;
    nCol = 0;
    nCount = 0;
    vectorNum = 0;
    inputVectorDegree = -1;
    degree = 1;
  }

  virtual ~PCA(){}

  Husky::ObjList<PCAObject>& getInputList() {
    return inputList;
  }

  int getDegree() {
    return degree;
  }

  void setDegree(int para){
    degree = para;
  }

  void loadFile(Husky::BaseWorker& worker, std::string fileNamePara) {
    Husky::HDFSLineInputFormat infmt;
    infmt.set_input(Husky::Context::get_params(fileNamePara));

    int marker = 0;

    worker.load(infmt, [&](boost::string_ref chunk) {
      if (chunk.empty()) return;
      boost::char_separator<char> sep(" \t");
      boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);

      for (auto &b : tok) {
        marker++;
        PCAObject tempPCAObject;
        if (marker % 2 == 0) {
          boost::char_separator<char> sep2(",");
          boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
       j   for (auto &c : tok2){
            tempPCAObject.inputTokens.push_back(c);
          }
          inputVectorDegree = tempPCAObject.size();
          worker.add_object(inputList, tempPCAObject);
        }
        else{
          tempPCAObject.label = std::stoi(b);
        }
      }
    });
  }

  void transform(){

    auto worker = Husky::Context::get_worker<Husky::Basewoker>();

    Husky::Aggregator<int> sum(0, [](int & a, const int & b){a += b});

    worker.list_execute(inputList, [&](PCAObject &tmpObj){
      sum.update(1);
    });

    vectorNum = sum.get_value();

    std::vector<double> initMean;
    initMean.resize(inputVectorDegree);//==========================斟酌一下

    Husky::Aggregator<std::vector<double>> meanVector(initMean, 
      [](std::vector<double> &a, const std::vector<double> &b) {
        for (int i = a.size() - 1; i >= 0; i--) a[i] += b[i];
      });

    worker.list_execute(inputList, [&] (PCAObject &tmpObj){
      meanVector.update(tmpObj.inputFeatures);
    });

    for (auto &x : meanVector.get_value()){
      x = x / vectorNum;
    }

    //step 1.零均值化

    worker.list_execute(inputList, [&] (PCAObject &tmpObj){
      for (int i = 0; i < inputVectorDegree; i++){
        tmpObj.inputFeatures[i] -= meanVector.get_value()[i];
      }
    });

    //step 2.求协方差矩阵
    //创造转置矩阵的Aggeragator

    std::vector<vector<double>> initTrans;
    initMean.resize(vectorNum);

    for (auto &x : initTrans){
      x.resize(inputVectorDegree);
    }

    Husky::Aggregator<std::vector<vector<double>>> transVector(initMean, 
      [](std::vector<vector<double>> &a, const std::vector<vector<double>> &b) {
        for (int i = a.size() - 1; i >= 0; i--)
          for (int j = i.size() - 1; j >= 0; j--){
            a[i][j] += b[i][j];
          } 
      },
      [=](std::vector<vector<double>> &a){
        a = initTrans;
      });

    //将值存入Aggeragator
    worker.list_execute(inputList, [&](PCAObject &tmpObj){
      int i = tmpObj.label;  // value at index 3 will be updated
      int len = tmpObj.inputFeatures.size();
      for (int j = 0; j < len; j++){
        transVector.update(std::make_pair(std::make_pair(i, j), val),
          [](std::vector<vector<double>> &a, const std::pair<std::pair<int, int>, double> & update) {
            std::pair<int, int> idx = update.first;
            double val = update.second;
            a[idx.second][idx.first] = val;  // you know these three lines can be written in one line, right
            });  // suppos
      }
    });

    //创造协方差矩阵Aggeragator
    /*std::vector<vector<double>> initCov;
    initCov.resize(inputVectorDegree);

    for (auto &x : initCov){
      x.resize(inputVectorDegree);
    }

    Husky::Aggregator<std::vector<vector<double>>> covMatrix(initCov, 
      [](std::vector<vector<double>> &a, const std::vector<vector<double>> &b) {
        for (int i = a.size() - 1; i >= 0; i--)
          for (int j = i.size() - 1; j >= 0; j--){
            a[i][j] += b[i][j];
          } 
      },
      [=](std::vector<vector<double>> &a){
        a = initCov;
      });*/

    //求协方差矩阵主程序
    for (int i = 0; i < transVector.get_value().size(); i++){
      covMatrixObject tempCovMatrixObject;
      tempCovMatrixObject.label = i;

      worker.list_execute(inputList, [&](PCAObject &tmpObj){
        int j = tmpObj.label;  // value at index 3 will be updated
        int len = tmpObj.inputFeatures.size();
        double val = 0;
        for (int k = 0; k < len; k++){
          val += transVector.get_value()[i][k] * tmpObj.inputFeatures[k];
        }

        tempCovMatrixObject.row[j] = val;
        /*covMatrix.update(std::make_pair(std::make_pair(i, j), val),
          [](std::vector<vector<double>> &a, const std::pair<std::pair<int, int>, double> & update) {
            std::pair<int, int> idx = update.first;
            double val = update.;
            a[idx.first][idx.secondsecond] = val;
            a[idx.second][idx.first] = val;  // you know these three lines can be written in one line, right
            });  // suppo*/
      });

      worker.add_object(covMatrix, tempCovMatrixObject);
    }

    //step 3.求特征值和特征向量

    /** 
* @brief 求实对称矩阵的特征值及特征向量的雅克比法  
* 利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量  
* @param pMatrix                长度为n*n的数组，存放实对称矩阵 
* @param nDim                   矩阵的阶数  
* @param pdblVects              长度为n*n的数组，返回特征向量(按列存储)  
* @param dbEps                  精度要求  
* @param nJt                    整型变量，控制最大迭代次数  
* @param pdbEigenValues         特征值数组 
* @return   
*/  
     //存储实对称矩阵的列向量，和行向量相同; 初始化特征向量矩阵
    worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj){
      for (auto &x : tmpObj.row){
        tmpObj.col.push_back(x);
      }

      int len = tmpObj.row.size();

      for (int i = 0; i < len; i++){
        if (i == tmpObj.label) {
          tmpObj.eigenvectorRow.push_back(1);
          tmpObj.eigenvectorCol.push_back(1);
        }
        else {
          tmpObj.eigenvectorRow.push_back(0);
          tmpObj.eigenvectorCol.push_back(0);
        }
      }
    });

    //开始迭代
    while(1) {
      dbMax = 0;
      worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj) {
        int len = tmpObj.eigenvectorRow.size();

        for (int i = 0; i < len; i++){
          if (fabs(tmpObj.eigenvectorRow[i]) > dbMax){
            dbMax = fabs(tmpObj.eigenvectorRow[i]);
            nRow = tmpObj.label;
            nCol = i;
          }
        }

        for (int i = 0; i < len; i++){
          if (fabs(tmpObj.eigenvectorCol[i]) > dbMax){
            dbMax = fabs(tmpObj.eigenvectorCol[i]);
            nCol = tmpObj.label;
            nRow = i;
          }
        }
      });

      if(dbMax < dbEps)     //精度符合要求   
        break;    

      if(nCount > nJt)       //迭代次数超过限制  
        break;  

      nCount++;

      worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj){
        if (tmpObj.label == nRow){
          dbApp = tmpObj.row[nRow];
          dbApq = tmpObj.row[nCol];
        }

        if (tmpObj.label = nCol){
          dbAqq = tmpObj.col[nCol];
        }
      });

      dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);  
      dbSinTheta = sin(dbAngle);  
      dbCosTheta = cos(dbAngle);  
      dbSin2Theta = sin(2*dbAngle);  
      dbCos2Theta = cos(2*dbAngle); 

      worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj){
        if (tmpObj.label == nRow){
          tmpObj.row[nRow] = dbApp*dbCosTheta*dbCosTheta + 
            dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;
          tmpObj.row[nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
        }

        if (tmpObj.label = nCol){
          tmpObj.col[nRow] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
          tmpObj.col[nCol] = dbApp*dbSinTheta*dbSinTheta +   
            dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;
        }
      });

      worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj){
        if ((tmpObj.label != nRow) && (tmpObj.label != nCol)){
          dbMax = tmpObj.row[nRow];
          tmpObj.row[nRow] = tmpObj.row[nCol]*dbSinTheta + dbMax*dbCosTheta; 
          tmpObj.row[nCol] = tmpObj.row[nCol]*dbCosTheta - dbMax*dbSinTheta;
        }
      });

      worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj){
        if ((tmpObj.label != nRow) && (tmpObj.label != nCol)){
          dbMax = tmpObj.col[nRow];
          tmpObj.col[nRow] = tmpObj.col[nCol]*dbSinTheta + dbMax*dbCosTheta;
          tmpObj.col[nCol] = tmpObj.col[nCol]*dbCosTheta - dbMax*dbSinTheta;
        }
      });

      //计算特征向量
      worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj){
          dbMax = tmpObj.eigenvectorRow[nRow];
          tmpObj.eigenvectorRow[nRow] = tmpObj.eigenvectorRow[nCol]*dbSinTheta + dbMax*dbCosTheta; 
          tmpObj.eigenvectorRow[nCol] = tmpObj.eigenvectorRow[nCol]*dbCosTheta - dbMax*dbSinTheta;
      });
    }

    //提取特征值
    worker.list_execute(covMatrix, [&](covMatrixObject &tmpObj){
      tmpObj.eigenValue = tmpObj.row[tmpObj.label];
    });
    //step 4.对特征值进行排序以及重新排列特征向量,特征值即pMatrix主对角线上的元素  
    auto topDegree = Husky::topk <double> (covMatrix,
        [](covMatrixObject& tmpObj){return tmpObj.eigenValue;},
        std::greater<double>(),
        degree);

    //output

    /*if (worker.id == 0) {
        Husky::log_msg("topk:");
        int rank = 1;
        for (const auto elem : topDegree) {
            Husky::log_msg(std::to_string(rank++)+"\t"+std::to_string(elem.first)+"\t"+elem.second);
        }
    }*/

    Husky::Context::free_worker<Husky::BaseWorker>();
  };
}

void runWorker() {
  Husky::Basewoker & worker = Husky::Context::get_worker<Husky::Basewoker>();
  auto & binarizerList = worker.create_list<BinarizerObject>("binarizer_list");

  Binarizer binarizer(binarizerList);

  binarizer.loadFile(worker, "input_binarizer");
  binarizer.transform(worker);
}

int main(int argc, char **argv) {
  Husky::run_job(runWorker, argv[1]);
  return 0;
}