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
using namespace std;
using namespace Husky;

class MatrixObject : public BaseObject{
public:
	int label;
	std::vector<double> row;

	typedef int KeyT;
	virtual KeyT const &id() const{
		return label;
	}
};

class BlockObject : public BaseObject{
public:
	pair<int, int> label;
	vector<vector<vecotr<double>>> block;

	typedef pair<int, int> KeyT;
	virtual KeyT const &id() const{
		return label;
	}
};

class MatrixCalculation {
private:
	int matrix1RowDegree;
	int matrix1ColDegree;
	int matrix2RowDegree;
	int matrix2ColDegree;
	int threadNum;
	int rowPerBlockForMatrix1;
	int colPerBlockForMatrix1;
	int rowPerBlockForMatrix2;
	int colPerBlockForMatrix2;
	int stepSize;
	int epochLength;

public:
	MatrixCalculation(int rowNum1, int colNum1, int rowNum2, int colNum2){
		matrix1ColDegree = colNum1;
		matrix1RowDegree = rowNum1;
		matrix2ColDegree = colNum2;
		matrix2RowDegree = rowNum2;
		/*threadNum = 0;
		getThreadNum();*/
		//calculate blocks for matrix1
		if (matrix1RowDegree % threadNum == 0){
			rowPerBlockForMatrix1 = matrix1RowDegree / rowNum1;
		}
		else if (matrix1RowDegree < threadNum){
			rowPerBlockForMatrix1 = matrix1RowDegree;
		}
		else{
			rowPerBlockForMatrix1 = matrix1RowDegree / colNum1+ 1;
		}

		if (matrix1ColDegree % threadNum == 0){
			colPerBlockForMatrix1 = matrix1ColDegree / rowNum2;
		}
		else if (matrix1ColDegree < threadNum){
			colPerBlockForMatrix1 = matrix1ColDegree;
		}
		else{
			colPerBlockForMatrix1 = matrix1ColDegree / colNum2 + 1;
		}

		//calculate blocks for matrix2
		if (matrix2RowDegree % threadNum == 0){
			rowPerBlockForMatrix2 = matrix2RowDegree / rowNum2;
		}
		else if (matrix2RowDegree < threadNum){
			rowPerBlockForMatrix2 = matrix2RowDegree;
		}
		else{
			rowPerBlockForMatrix2 = matrix2RowDegree / rowNum2 + 1;
		}

		if (matrix2ColDegree % threadNum == 0){
			colPerBlockForMatrix2 = matrix2ColDegree / rowNum2;
		}
		else if (matrix2ColDegree < threadNum){
			colPerBlockForMatrix2 = matrix2ColDegree;
		}
		else{
			colPerBlockForMatrix2 = matrix2ColDegree / colNum2 + 1;
		}
	}

	virtual ~MatrixCalculation(){}

	//load matrix1 from file
	void loadMatrix1(ObjList &matrix, string fileNamePara) {
		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(fileNamePara));

    	worker.load(infmt, [&](boost::string_ref chunk) {
      		if (chunk.empty()) return;
      		boost::char_separator<char> sep(" \t");
      		boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);

      		for (auto &b : tok) {
        		MatrixObject tempMatrixObj;
          		boost::char_separator<char> sep2(",");
          		boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);

          		for (auto &c : tok2){
          			tempMatrixObj.row.push_back(stod(c));
          		}
          	}

          	matrix1RowDegree = tempPCAObject.size();
          	tempMatrixObj.label = worker.get_size(matrix) + 1;
          	worker.add_object(matrix, tempMatrixObj);
        });
        matrix1ColDegree = worker.get_size(matrix) + 1;

	}

	//load matirx2 from file
	void loadMatrix2(ObjList &matrix, string fileNamePara) {
		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(fileNamePara));

    	worker.load(infmt, [&](boost::string_ref chunk) {
      		if (chunk.empty()) return;
      		boost::char_separator<char> sep(" \t");
      		boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);

      		for (auto &b : tok) {
        		MatrixObject tempMatrixObj;
          		boost::char_separator<char> sep2(",");
          		boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);

          		for (auto &c : tok2){
          			tempMatrixObj.row.push_back(stod(c));
          		}
          	}

          	matrix2RowDegree = tempPCAObject.size();
          	tempMatrixObj.label = worker.get_size(matrix) + 1;
          	worker.add_object(matrix, tempMatrixObj);
        });
        matrix2ColDegree = worker.get_size(matrix) + 1;

	}

	void getThreadNum(){
		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params("socket_file"));

    	worker.load(infmt, [&](boost::string_ref chunk) {
      		if (chunk.empty()) return;
      		boost::char_separator<char> sep(" \t");
      		boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);

      		for (auto &b : tok) {
        		MatrixObject tempMatrixObj;
          		boost::char_separator<char> sep2(":");
          		boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
          		auto it = tok2.begin();
          		it++;
          		threadNum += stoi(*it++);
          	}
        });
	}

	void RowMatrix2BlockMatrix(ObjList &rowMatrix, ObjList &blockMatrix, int rowPerBlock, int colPerBlock){
		worker.list_execute(rowMatrix, [&](MatrixObj &tmpObj) {
			int len = tmpObj.row.size();
			vector<double> tempPartitonVec;
			for (int i = 0; i < len; i++){
				if (i % colPerBlock == 0){
					tempPartitonVec.clear();
					tempPartitonVec.push_back(i / colPerBlock);
				}
				tempPartitonVec.push_back(tmpObj.row[i]);
				if (i % colPerBlock == colPerBlock - 1){
					 worker.send_message<vector<double>>(tempPartitonVec, make_pair(tmpObj.label / rowPerBlock, i / colPerBlock), w_list);
				}
			}
		});

		worker.register_msg_ctor<vector<double>>(blockMatrix, [] (vector<double> msg, BlockObject::KeyT key){
			BlockObject n; 
			n.label = key; 
			return n; 
		});

    	worker.list_execute(blockMatrix, [&](BlockObject &tmpObj) {
    		auto tempBlock = tmpObj.get_messages<vector<double>>();
    		tmpObj.block.resize(tempBlock.size());
    		for (auto &x : tempBlock){
    			int len = x.size();
    			tmpObj.block[x[0]].resize(len - 1);
    			for (int k = 1; k < len; k++){
    				tmpObj.block[x[0]][k - 1] = x[k];
    			}
    		}
    	});

    	//free row matrix
    	worker.free_list(rowMatrix);
	}

	void mulBlock(BlockObject &block1, BlockObject &block2, BlockObject &res){

		int len1 = block1.block.size();
		int len3 = block2.block[0].size();
		int len2 = block2.block.size();

		res.block.resize(len1);
		for (auto &x : res.block){
			x.resize(len3);
		}

		for (int i = 0; i < len1; i++)
			for (int j = 0; j < len3; j++){
				res.block[i][j]=0;
				for (int k = 0; k < len2; k++)
					res.block[i][j] += block1.block[i][k] * block2.block[k][j];             
			}    
	}

	//spark要求最少能存下至少左矩阵一行和右矩阵一列的非零矩阵块，所以可以将消息发送给
	void mulBlockMatrix(ObjList &blockList1, ObjList &blockList2, ObjList &res){

		worker.list_execute(blockList1, [&](BlockObject &block1) {

				worker.request(doc.words.at(i), term_list);
		});
    worker.list_reply(term_list, [](Term& t){
        return t.idf;
    });

    worker.list_execute(document_list, [&](Document& doc) {
        for (int i = 0; i < doc.tf_idf.size(); i++)
            doc.tf_idf.at(i) = doc.tf_idf.at(i)*doc.get_response<double>(doc.words.at(i));
    });
	}

	void add(ObjList &matrix, double num) {
		worker.list_execute(matrix, [&](MatrixObject &tmpObj)){
			for (auto &x : tmpObj.row) {
				x += num;
			}
		}
	}

	void add(ObjList &matrix, int row, int col, double num){
		MatrixObject &x = matrix.find(row);
		if (x != nullptr){
			if (col <= x.row.size()){
				x.row[col] += num;
			}
			else{
				log_msg("Wrong Column!");
			}
		}
		else{
			log_msg("Wrong Row!");
		}
	}

	void add(ObjList &matrix1,  ObjList &matrix2){
		worker.list_execute(matrix1, [&](MatrixObject &tmpObj) {
			worker.request(tmpObj.id(), matrix2);
		});

		worker.list_reply(matrix2, [](MatrixObject &tmpObj){
			return tmpObj.row;
		});

		worker.list_execute(matrix1, [&](MatrixObject &tmpObj) {
			auto vec = tmpObj.get_response<vector<double>>(tmpObj.id());
			int len = tmpObj.row.size();

			for (int i = 0; i < len; i++){
				tmpObj.row[i] += vec[i];
			}
		});
	}

	void sub(ObjList &matrix, double num) {
		worker.list_execute(matrix, [&](MatrixObject &tmpObj)){
			for (auto &x : tmpObj.row) {
				x -= num;
			}
		}
	}

	void sub(ObjList &matrix, int row, int col, double num){
		MatrixObject &x = matrix.find(row);
		if (x != nullptr){
			if (col <= x.row.size()){
				x.row[col] -= num;
			}
			else{
				log_msg("Wrong Column!");
			}
		}
		else{
			log_msg("Wrong Row!");
		}
	}

	void sub(ObjList &matrix1,  ObjList &matrix2){
		worker.list_execute(matrix1, [&](MatrixObject &tmpObj) {
			worker.request(tmpObj.id(), matrix2);
		});

		worker.list_reply(matrix2, [](MatrixObject &tmpObj){
			return tmpObj.row;
		});

		worker.list_execute(matrix1, [&](MatrixObject &tmpObj) {
			auto vec = tmpObj.get_response<vector<double>>(tmpObj.id());
			int len = tmpObj.row.size();

			for (int i = 0; i < len; i++){
				tmpObj.row[i] -= vec[i];
			}
		});
	}

	void mul(ObjList &matrix, double num){
		worker.list_execute(matrix, [&](MatrixObject &tmpObj)){
			for (auto &x : tmpObj.row) {
				x *= num;
			}
		}
	}

	void mul(ObjList &matrix, int row, int col, double num){
		MatrixObject &x = matrix.find(row);
		if (x != nullptr){
			if (col <= x.row.size()){
				x.row[col] -= num;
			}
			else{
				log_msg("Wrong Column!");
			}
		}
		else{
			log_msg("Wrong Row!");
		}
	}

	void mul(ObjList &matrix1,  ObjList &matrix2, ObjList &blockList1, ObjList &blockList2, ObjList &res/*, ObjList &partitionList1, ObjList &partitionList2*/){

		//when error
		if (matrix1ColDegree != matrix2RowDegree){
			log_msg("These two matricies can't be multiplied");
			return;
		}

		//split whole matrix into blocks
		RowMatrix2BlockMatrix(matrix1, blockList1, rowPerBlockForMatrix1, colPerBlockForMatrix1);
		RowMatrix2BlockMatrix(matrix2, blockList2, rowPerBlockForMatrix2, colPerBlockForMatrix2);

		//mul blocks
		mulBlockMatrix(blockList1, blockList2, res);
	}

	void pow(ObjList &matrix, ObjList &res, int num){
		res = matrix;
		ObjList<BlockObject> &blockList1;
		ObjList<BlockObject> &blockList2;
		for (int i = 1; i < num; i++){
			mul(res, matrix, blockList1, blockList2, res);
		}
	}

	void transpose(ObjList &matrix,  ObjList &resMatrix){
		worker.list_execute(matrix, [&](MatrixObj &tmpObj) {
			int len = tmpObj.row.size();
			for (int i = 0; i < len; i++){
				worker.send_message<pair<int, double>>(make_pair(tmpObj.label, tmpObj.row[i]), i, resMatrix);
			}
		});

		worker.register_msg_ctor<pair<int, double>>(resMatrix, [] (pair<int, double> msg, MatrixObject::KeyT key){
			BlockObject n; 
			n.label = key; 
			return n; 
		});

    	worker.list_execute(resMatrix, [&](MatrixObject &tmpObj) {
    		auto tempRow = tmpObj.get_messages<pair<int, double>>();
    		tmpObj.row.resize(tempRow.size());
    		for (auto &x : tempRow){
    			tmpObj.row[x.first] = x.second;
    		}
    	});	
	}

	//approximate
	void engineVectorApproximate(ObjList &matrix,  ObjList &engineVectorMatrix, vector<double> &resVector){
		//1.先求出PCA的最开始的最重要的向量
		//创造转置矩阵
		auto & transposeMatrix = worker.create_list<MatrixObject>("transposeMatrix");
		transpose(matrix, transposeMatrix);
		//统计整个矩阵的列数
		Husky::Aggregator<int> rowSum(0, [](int & a, const int & b){ a += b; });
		worker.list_execute(matrix, [&](MatrixObject &tmpObj){
			rowSum.update(1);
		});
		//统计整个矩阵的行数
		Husky::Aggregator<int> colSum(0, [](int & a, const int & b){ a += b; });
		worker.list_execute(transposeMatrix, [&](MatrixObject &tmpObj){
			colSum.update(1);
		});


		//创造波浪U的Aggregator
		std::vector<double> initTidleU;
		tidleU.resize(rowSum.get_value());//列向量的长度

		Husky::Aggregator<std::vector<double>> tidleU(initTidleU,
			[](std::vector<double> &a, const std::vector<double> &b){
				int len = a.size();
				for (int i = 0; i < len; i++){
					a[i] += b[i];
				}
			});
		std::vector<double> v;
		//init tidle w

		//for loop
		for (int s = 0; s < stepSize; s++){
			//1.更新波浪u
			tidleU.value = initTidleU;
			worker.list_execute(transposeMatrix, [&](MatrixObject &tmpObj){
				double mul = 0;
				int len = tmpObj.row.size();

				for (int i = 0; i < len; i++){
					mul += tmpObj.row[i] * tidleW[i];
				}

				std::vector<double> tempRow;
				for (int i = 0; i < len; i++){
					tempRow.push_back(x[i] * mul / colSum.get_value());
				}

				tidleU.update(tempRow);
			});

			auto w0 = tidleW.get_value();
			std::vector<double> wt;

			for (int t = 0; t < epochLength; t++){
				int k = rand();
				auto Xit = transposeMatrix.find(k).row;
				int len = xk.size();
				double mul1 = 0, mul2 = 0;
				auto tidleWtSub1 = tidleW.get_value();
				for (int i = 0; i < len; i++){
					mul1 += Xit[i] * w0[i];
					mul2 += Xit[i] * tidleWtSub1[i];
				}

				std::vector<double> tempWt;
				double dis = 0;
				for (int i = 0; i < len; i++){
					double a = wt[i] + stepSize * (Xit[i] * (mul1 - mul2) + tidleU.get_value()[i]);
					tempWt[i] = a;
					dis += a * a;
				}
				dis = sqrt(dis);

				for (int i = 0; i < len; i++){
					wt[i] = (1 / dis) * tempWt[i];
				}
			}

			for (int i = 0; i < tidleW.size){
				tidleW[i] = wt[i];
				v.push_back(wt[i]);
			}			
		}

		//2.根据第一个最重要的向量推导以后的
		//创造1/n乘一大堆x
		std::vector<double> initDeflatedXMatrix;
		initDeflatedXMatrix.resize(rowSum.get_value());//列向量的长度

		Husky::Aggregator<std::vector<double>> deflatedXMatrix(initDeflatedXMatrix,
			[](std::vector<double> &a, const std::vector<double> &b){
				int len = a.size();
				for (int i = 0; i < len; i++){
					a[i] += b[i];
				}
			});
		worker.list_execute(transposeMatrix, [&](MatrixObject &tmpObj){
				double mul = 0;
				int len = tmpObj.row.size();

				std::vector<double> tempRow;
				for (int i = 0; i < len; i++){
					tempRow.push_back(x[i] * x[i] / colSum.get_value());
				}

				deflatedXMatrix.update(tempRow);
			});

		int pos = 0;

		//output engineVector1
		int vSize = v.size();
		std::string log = std::to_string(pos) + " : [";
		for (int i = 0; i < vSize; i++){
			log += std::to_string(v[i]);
			log += ",";
		}
		log += "]";
		Husky::log_msg(log);

		std::vector<double> vSum;
		for (int i = 0; i < vSize; i++){
				vSum.push_back(0);
		}

		//calculate other engineVectors
		while (pos < colSum.get_value()){
			pos++;
			for (int i = 0; i < vSize; i++){
				vSum[i] += v[i];
			}

			//calculate and output
			auto tempX = deflatedXMatrix.get_value();
			std::string log = std::to_string(pos) + " : [";
			for (int i = 0; i < vSize; i++){
				v[i] = tempX[i] - vSum[i];
				log += std::to_string(v[i]);
				log += ",";
			}
			log += "]";
			Husky::log_msg(log);
		}
	}

	void engineVectorAndEngineValueExact(int nDim, double *pdblVects, double *pdbEigenValues, double dbEps, int nJt){
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
		double * pMatrix;

		worker.list_execute(matrix, [&](MatrixObject &tmpObj){
			int len = tmpObj.row.size();
			for (int i = 0; i < len; i++){
				pMatrix[tmpObj.label * nDim + i] = tmpObj.row[i];
			}
		});

		for(int i = 0; i < nDim; i ++)   
		{     
			pdblVects[i*nDim+i] = 1.0f;   
			for(int j = 0; j < nDim; j ++)   
			{   
				if(i != j)     
					pdblVects[i*nDim+j]=0.0f;   
			}   
		}   

		int nCount = 0;     //迭代次数  
		while(1)  
		{  
			//在pMatrix的非对角线上找到最大元素  
			double dbMax = pMatrix[1];  
			int nRow = 0;  
			int nCol = 1;  
			for (int i = 0; i < nDim; i ++)          //行  
			{  
				for (int j = 0; j < nDim; j ++)      //列  
				{  
					double d = fabs(pMatrix[i*nDim+j]);   
					if((i!=j) && (d> dbMax))   
					{   
						dbMax = d;     
						nRow = i;     
						nCol = j;   
					}   
				}  
			}  

			if(dbMax < dbEps)     //精度符合要求   
            	break;    
  
        	if(nCount > nJt)       //迭代次数超过限制  
            	break;  
  
        	nCount++;  
  
        	double dbApp = pMatrix[nRow*nDim+nRow];  
        	double dbApq = pMatrix[nRow*nDim+nCol];  
        	double dbAqq = pMatrix[nCol*nDim+nCol];  
  
        	//计算旋转角度  
        	double dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);  
        	double dbSinTheta = sin(dbAngle);  
        	double dbCosTheta = cos(dbAngle);  
        	double dbSin2Theta = sin(2*dbAngle);  
        	double dbCos2Theta = cos(2*dbAngle);  
  
        	pMatrix[nRow*nDim+nRow] = dbApp*dbCosTheta*dbCosTheta +   
            dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;  
        	pMatrix[nCol*nDim+nCol] = dbApp*dbSinTheta*dbSinTheta +   
            dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;  
        	pMatrix[nRow*nDim+nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;  
        	pMatrix[nCol*nDim+nRow] = pMatrix[nRow*nDim+nCol];  
  
        	for(int i = 0; i < nDim; i ++)   
        	{   
            	if((i!=nCol) && (i!=nRow))   
            	{   
                	int u = i*nDim + nRow;  //p    
                	int w = i*nDim + nCol;  //q  
                	dbMax = pMatrix[u];   
                	pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta;   
                	pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta;   
            	}   
        	}   
  
        	for (int j = 0; j < nDim; j ++)  
        	{  
            	if((j!=nCol) && (j!=nRow))   
            	{   
                	int u = nRow*nDim + j;  //p  
                	int w = nCol*nDim + j;  //q  
                	dbMax = pMatrix[u];   
                	pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta;   
                	pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta;   
            	}   
        	}  
  
        	//计算特征向量  
        	for(int i = 0; i < nDim; i ++)   
        	{   
            	int u = i*nDim + nRow;      //p     
            	int w = i*nDim + nCol;      //q  
            	dbMax = pdblVects[u];   
            	pdblVects[u] = pdblVects[w]*dbSinTheta + dbMax*dbCosTheta;   
            	pdblVects[w] = pdblVects[w]*dbCosTheta - dbMax*dbSinTheta;   
        	}   
  
    	}  
  
    	//对特征值进行排序以及重新排列特征向量,特征值即pMatrix主对角线上的元素  
    	std::map<double,int> mapEigen;  
    	for(int i = 0; i < nDim; i ++)   
    	{     
        	pdbEigenValues[i] = pMatrix[i*nDim+i];  
  
        	mapEigen.insert(make_pair( pdbEigenValues[i],i ) );  
    	}   
  
    	double *pdbTmpVec = new double[nDim*nDim];  
    	std::map<double,int>::reverse_iterator iter = mapEigen.rbegin();  
    	for (int j = 0; iter != mapEigen.rend(),j < nDim; ++iter,++j)  
    	{  
        	for (int i = 0; i < nDim; i ++)  
        	{  
            	pdbTmpVec[i*nDim+j] = pdblVects[i*nDim + iter->second];  
        	}  
  
        	//特征值重新排列  
        	pdbEigenValues[j] = iter->first;  
    	}  
  
    	//设定正负号  
    	for(int i = 0; i < nDim; i ++)   
    	{  
        	double dSumVec = 0;  
        	for(int j = 0; j < nDim; j ++)  
            	dSumVec += pdbTmpVec[j * nDim + i];  
        	if(dSumVec<0)  
        	{  
            	for(int j = 0;j < nDim; j ++)  
                	pdbTmpVec[j * nDim + i] *= -1;  
        	}  
    	}  
  
    	memcpy(pdblVects,pdbTmpVec,sizeof(double)*nDim*nDim);  
    	delete []pdbTmpVec;    
	}
};

void runWorker() {
	Husky::Basewoker & worker = Husky::Context::get_worker<Husky::Basewoker>();
	auto & matrix1 = worker.create_list<MatrixObject>("matrix1");
	auto & matrix2 = worker.create_list<MatrixObject>("matrix2");
	/*auto & partitionList1 = worker.create_list<MatrixObject>("partition_list1");
	auto & partitionList2 = worker.create_list<MatrixObject>("partition_list2");*/
	auto & blockList = worker.create_list<BlockObject>("block_list");
}

int main(int argc, char **argv) {
	Husky::run_job(runWorker, argv[1]);
	return 0;
}