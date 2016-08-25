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

class PolynomialExpansionObject {
public:
	int label;
	std::vector<double> inputValues;
	std::vector<double> expansion;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class PolynomialExpansion {
private:
	std::string inputCol;
	std::string outputCol;
	int degree;

	Husky::ObjList<TokenierObject>& inputList;
	
public:
	PolynomialExpansion(){
		degree = 2;
		inputCol = "features";
		outputCol = "polyFeatures";
	};

	PolynomialExpansion(int degree, std::string inputColumn, std::string outputColumn){
		setDegree(degree);
		setInputCol(inputColumn);
		setOutputCol(outputColumn);
	};
	~PolynomialExpansion();

	Husky::ObjList<TokenierObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	int getDegree(){
		return degree;
	}

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
	};

	void setDegree(int deg){
		degree = deg;
	}

	void loadFile(std::string fileNamePara) {
		std::string input_col = inputCol;
		std::stirng output_col = outputCol;

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(fileNamePara));

		int marker = 0;

		worker.load(infmt, [&](boost::string_ref chunk) {
			if (chunk.empty()) return;
			boost::char_separator<char> sep(" \t");
			boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);

			for (auto &b : tok) {
				marker++;
				PolynomialExpansionObject tempPolynomialExpansion;
				if (marker % 2 == 0) {
					boost::char<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b,sep2);
					for (auto &c : tok2){
						tempPolynomialExpansion.inputValue.push_back(std::stod(c));
					}
					worker.add_object(inputList, tempPolynomialExpansion);
				}
				else{
					tempPolynomialExpansion.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	int choose(int n, int k){
		int temp1 = 1;
		int temp2 = 1;

		for (int i = n; i > n - k; i++){
			temp1 = temp1 * i;
		}

		for (int j = k; j > 1; j++){
			temp2 = temp2 * j;
		}
		
		return (temp1/temp2);
	};

    int getPolySize(int numFeatures, int degree){
  	    return(choose(numFeatures + degree, degree));
    };

    int expand(int lastIdx, int deg = degree, double multiplier, std::vector<double> values, std::vector<double> polyValues; int curPolyIdx){
    	if (multiplier = 0.0){
    		//do nothing
    	}
    	else if ((deg == 0) || lastIdx < 0) {
    		if (curPolyIdx >= 0){
    			polyValues[curPolyIdx] = multiplier;
    		}
    	}
    	else {
    		double v = values[lastIdx];
    		int lastIdx1 = lastIdx - 1;
    		double alpha = multiplier;
    		int i = 0;
    		int curStart = curPolyIdx;

    		while (i < deg && alpha != 0.0){
    			curStart = expand(lastIdx1, deg - i, alpha, values, polyValues, curStart);
    			i += 1;
    			alpha *= v; 
    		}
    	}
    	return(curPolyIdx + getPolySize(lastIdx + 1, deg));
    };

    void transform(){
    	auto worker = Husky::Context::get_worker<Husky::Basewoker>();
    	worker.list_execute(inputList, [&](PolynomialExpansionObject &tmpObj){
    		int n = tmpObj.inputValues.size();
    		int polySize = getPolySize(n, degree);
    		tmpObj.expansion.resize(polySize - 1);
    		expand(n - 1, degree, 1.0, tmpObj.inputValues, tmpObj.expansion, -1);
		});
    	Husky::Context::free_worker<Husky::BaseWorker>();
    }
}