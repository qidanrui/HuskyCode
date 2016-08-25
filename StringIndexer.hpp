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

class StringIndexerObject{
public:
	int label;
	std::string category;
	double index;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}

};

class WordCountObject{
public:
	std::string label;
	int count;
	typedef std::string KeyT;

	virtual KeyT const & id() const {
		return label;
	}	
};

class StringIndexer {
private:
	std::string inputCol;
	std::string outputCol;

	Husky::ObjList<StringIndexerObject>& inputList;
	Husky::ObjList<WordCountObject> & wordCountList;
	
public:
	StringIndexer(){
		inputCol = "inputTokens";
		outputCol = "stringIndexer";
	};

	StringIndexer(int para, std::string inputColumn, std::string outputColumn){
		inputCol = inputColumn;
		outputCol = outputColumn;
	};
	~StringIndexer();

	Husky::ObjList<StringIndexerObject>& getInputList() {
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
				StringIndexerObject tempStringIndexerObject;
				if (marker % 2 == 0) {
					tempStringIndexerObject.category = b;
					worker.add_object(inputList, tempStringIndexerObject);
					worker.send_message<Husky::SumCombine<int>>(1, b, wordCountList);
				}
				else{
					tempNGramObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();	

		Husky::Aggregator<int> sum(0, [](int &a, const int & b){
			a += b;
		});

		worker.register_msg_ctor<int>(wordCountList, [](int msg, WordCountObject::KeyT key){
			WordCountObject n;
			n.label = key;
			return n;
		});

		worker.list_execute(wordCountList, [&](WordCountObject & w){
			w.count = std::accumulate(w.get_messages<int>().begin(), w.get_messages<int>().end(), 0);
			sum.update(1);
		});


		/*std::unordered_map<std::string, std::pair(int, double)> initWordLabel;
		Husky::Aggregator<std::unordered_map<std::string, std::pair(int, double)>> wordLabel(initWordLabel,
			[](std::unordered_map<std::string, std::pair(int, double)> &a, const std::unordered_map<std::string, std::pair(int, double)> &b){
				for (int i = 0; i < )				
			},
			[&](std::unordered_map<std::string, double> &a){
				a = initWordLabel;
			});
		worker.list_execute(inputList, [&](StringIndexerObject &tmpObj){
			std::string index = tmpObj.catagory;
			double val = tmpObj.catagory;
			unordered_map.update(std::make_pair(index, val),
				[](std::unordered_map<std::string &a, const std::pair<std::string, double> &update){
					std::string idx = update.first;
					double val = update.second;
					a[idx] += val;
				});
		});*/

		auto sortResult = Husky::topk <int> (w_list, [](WordCountObject& w){return w.num_occur;}, std::greater<int>(), sum.get_value());
		//对Aggeragator排序

		worker.list_execute(inputList, [&](StringIndexerObject &tmpObj){
			tmpObj.index = double(sortResult[tmpObj.catagory]);
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}