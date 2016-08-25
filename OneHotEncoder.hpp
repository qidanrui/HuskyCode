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

class OneHotEncoderObject{
public:
	int label;
	std::vector<double> inputValues;
	std::vector<double> encoder;
	typedef int KeyT;

	virtual KeyT const & id() const{
		return label;
	}
};

class CatagoryObject{
public:
	int label;
	std::vector<std::pair<std::string, int>> categoryCount;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
};

class  OneHotEncoder{
private:
	std::string inputCol;
	std::string outputCol;

	Husky::ObjList<OneHotEncoderObject>& inputList;
	Husky::ObjList<CatagoryObject> & catagoryList;
	
public:
	OneHotEncoder(){
		inputCol = "inputTokens";
		outputCol = "stringIndexer";
	};

	OneHotEncoder(std::string inputColumn, std::string outputColumn){
		inputCol = inputColumn;
		outputCol = outputColumn;
	};
	~OneHotEncoder();

	Husky::ObjList<OneHotEncoderObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
	};

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
				OneHotEncoderObject tempOneHotEncoderObject;

				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempOneHotEncoderObject.inputValues.push_back(c);
						worker.send_message<std::string>(c, &c - &tok2, catagoryList);
					}
					worker.add_object(inputList, tempOneHotEncoderObject);
				}
				else{
					tempOneHotEncoderObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};


    bool SortByM1( const std::pair<std::string, int> &v1, const std::pair<std::string, int> &v2){  
    	return v1.second > v2.second;//升序排列  
    };

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.register_msg_ctor<int>(wordCountList, [](int msg, WordCountObject::KeyT key){
			WordCountObject n;
			n.label = key;
			return n;
		});

		worker.list_execute(CatagoryList, [&](CatagoryObject & w){
			for (auto &x : w.get_messages<std::string>()){
				for (auto &b : w.catagoryCount) {
					if (x == b.first) b.second++;
					else {
						w.categoryCount.push_back(std::make_pair(b.first, 1));
					}
				}
			}
			std::sort(w.categoryCount.begin(),w.categoryCount.end(),SortByM1);  
		});

		worker.globalize_list(CatagoryList);

		worker.list_execute(inputList, [&](OneHotEncoderObject &tmpObj){
			len = tmpObj.inputValues.size();

			for (int i = 0; i < len; i++){
				for (int j = 0; j < catagoryList[i].catagoryCount.size(); j++){
					if (tmpObj.inputValues[i] == catagoryList[i].categoryCount[j].first){
						for (int k = 0; k < catagoryList[i].categoryCount[j].first - 1; k++){
							tmpObj.encoder.push_back(0);
						}
						tmpObj.encoder[-1] = 1;
						for (int k = catagoryList[i].categoryCount[j].first - 1; k <catagoryList[i].categoryCount.size(); k++){
							tmpObj.encoder.push_back(0);
						}
					}
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}