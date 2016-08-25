// Copyright 2015 Husky Team
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <functional>
#include <iostream>
#include <string>

#include "boost/tokenizer.hpp"
#include "gperftools/profiler.h"

#include "core/engine.hpp"
#include "io/inputformat/lineinputformat.hpp"
#include "lib/topk.hpp"


class vocabWordObject : public Husky::BaseObject {
public:
    std::string word;
    long long num_occur;//词频
    std::string code; //huffman编码
    int *point;//huffman编码对应内节点的路径  
    int codelen;
    typedef std::string KeyT;

    virtual KeyT const & id() const {
        return word;
    }
};

//对输入的train_data或者output_data做wordCount
void wordCount(Huskey::ObjList<vocabWordObject> *w_list) {
    Husky::BaseWorker & worker = Husky::Context::get_worker<Husky::BaseWorker>();

    auto & w_list = worker.create_list<vocabWordObject>("words");

    auto parse_wc = [&] (boost::string_ref & chunk) {
        if (chunk.size() == 0) return;
        boost::char_separator<char> sep(" \t");
        boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);
        for (auto &w : tok) {
            worker.send_message<Husky::SumCombine<int>>(1, w, w_list);
        }
    };

    #ifdef WITH_HDFS8
        Husky::LineInputFormat<Husky::HDFSFileSplitter> infmt;
    #else
        Husky::LineInputFormat<Husky::LocalFileSplitter> infmt;
    #endif
    infmt.set_input(Husky::Context::get_params("input"));
    worker.load(infmt, parse_wc);

    worker.register_msg_ctor<int>(w_list, [] (int msg, vocabWordObject::KeyT key){ vocabWordObject n; n.word = key; return n; });

    worker.list_execute(w_list, [&](vocabWordObject & w) {
        w.num_occur = std::accumulate(w.get_messages<int>().begin(), w.get_messages<int>().end(), 0);
    });

    Husky::Context::free_worker<Husky::BaseWorker>();
}

//搜索一个词并且返回一个指针
vocabWordObject & SearchVocab(Husky::ObjList<vocabWordObject> &vocabList, std::string word)  
{  
  return(vocabList.find(word));
}  

//从本地读取词汇表
void ReadVocab(Husky::ObjList<vocabWordObject> &vocabList)  
{  
  Husky::BaseWorker & worker = Husky::Context::get_worker<Husky::BaseWorker>();

  auto parse_wc = [&] (boost::string_ref & chunk) {
      if (chunk.size() == 0) return;
      boost::char_separator<char> sep(" \t");
      boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);
      for (auto &w : tok) {
          boost::char_separator<char> sep2(':');
          boost::tokenizer<boost::char_separator<char>> tok2(chunk, sep2);
          auto it = tok2.begin();
          vocabWordObject tempVocabWordObject;
          tempVocabWordObject.word = *it;
          tempVocabWordObject.num_occur = *it++;
          worker.add_object(vocabList, tempVocabWordObject);
      }
  };
  #ifdef WITH_HDFS8
      Husky::LineInputFormat<Husky::HDFSFileSplitter> infmt;
  #else
      Husky::LineInputFormat<Husky::LocalFileSplitter> infmt;
  #endif
  infmt.set_input(Husky::Context::get_params("inputVocabList"));
  worker.load(infmt, parse_wc);
}  

//将train data读取到的词加到词汇表中
void AddWordToVocab(std::string word)  
{  
  //proprecssion
  unsigned int hash, length = strlen(word) + 1;  
  if (length > MAX_STRING)  
      length = MAX_STRING;  
  //add word to vocabulary
  vocabWordObject tempVocabWordObject;
  tempVocabWordObject.word = word;
  tempVocabWordObject.num_occur = 1;
}  

//若词汇表太大则缩减
void ReduceVocab(Husky::ObjList<vocabList> &vocabList){
  worker.list_execute(vocabList, [&](vocabWordObject &tmpObj){
    if (tmpObj.num_occur < min_reduce){
      worker.delete_object(vocabList, tmpObj);
    }
  });
}

//根据词频排序
void SortVocab(Husky::ObjList<vocabWordObject> &vocabList){
  //先移除词频过小的词
  worker.list_execute(vocabList, [&](vocabWordObject &tmpObj){
    if (tmpObj.num_occur < min_count){
      worker.delete_object(vocabList, tmpObj);
    }
  });

  //进行排序===================
}

//从Train data中读取词汇
void LearnVocabFromTrainFile(Husky::ObjList<vocabWordObject> &vocabList)  
{  
  wordCount(vocabList);
  SortVocab(vocabList);
}  

//对网络进行初始化并且创造哈夫曼树
void InitNet()  
{  
  /*long long a, b;  
  a = posix_memalign((void **)&syn0, 128, (long long)vocab_size * layer1_size * sizeof(real));  
  //先知道这个也是申请动态数组，对齐还有128这个参数以后再了解  
  if (syn0 == NULL)  
  {  
      printf("Memory allocation failed\n"); exit(1);  
  }  
  if (hs)//采用softmax  
  {  
    a = posix_memalign((void **)&syn1, 128, (long long)vocab_size * layer1_size * sizeof(real));  
    if (syn1 == NULL)  
    {  
        printf("Memory allocation failed\n"); exit(1);  
    }  
    for (b = 0; b < layer1_size; b++)  
        for (a = 0; a < vocab_size; a++)  
            syn1[a * layer1_size + b] = 0;  
  }  
  if (negative>0)//还有负样本  
  {  
    a = posix_memalign((void **)&syn1neg, 128, (long long)vocab_size * layer1_size * sizeof(real));  
    if (syn1neg == NULL)  
    {  
        printf("Memory allocation failed\n"); exit(1);  
    }  
    for (b = 0; b < layer1_size; b++)  
        for (a = 0; a < vocab_size; a++)  
            syn1neg[a * layer1_size + b] = 0;  
  }  
  for (b = 0; b < layer1_size; b++)  
      for (a = 0; a < vocab_size; a++)  
          syn0[a * layer1_size + b] = (rand() / (real)RAND_MAX - 0.5) / layer1_size;  */
  CreateBinaryTree();//建立huffman树，对每个单词进行编码  
}  

bool SortByM1( const pair<std::string, int> &v1, const pair<std::string, int> &v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
{  
    return v1.second < v2.second;//升序排列  
}  

void CreateBinaryTree (Husky::ObjList<wordObject> &vocabList) {
  //得到词库的词的总数
  Husky::Aggregator<int> sum(0, [](int & a, const int & b){ a += b; });
  worker.list_execute(vocabList, [&](vocabWordObject &tmpObj){
    sum.update(1);
  }]);

  //以Aggregator创造哈夫曼树
  std::vector<pair<std::string, int>> initHuffmanTree;
  initHuffmanTree.resize(2 * sum.get_value() - 1);

  Husky::Aggregator<std::vector<pair<std::string, int>>> huffmanTree(initHuffmanTree,
    [](std::vector<pair<std::string, int>> &a, const std::vector<pair<std::stringm int>> &b) {
      a[a.size() + 1] = b[0];
    });

  worker.list_execute(vocabList, [&](vocabWordObject &tmpObj){
    huffmanTree.update(std::make_pair(tmpObj.word, tmpObj.num_occur));
  });

  /*worker.list_execute(vocabList, [&](vocabWordObject &tmpObj){

  });*/

  if (worker.id == 0){
    auto tempHuffmanTree = huffmanTree.get_value();
    int len = sum.get_value();
    int mark = len;
    int i;

    long long *binary = (long long *)calloc(len * 2 + 1, sizeof(long long));  
    long long *parent_node = (long long *)calloc(len * 2 + 1, sizeof(long long));  
    

    int pos1 = len - 1;
    int pos2 = len;

    int min1i, min2i;

    for (i = 0; i < len; i++)  
    {  
        // First, find two smallest nodes 'min1, min2' 找出目前权值最小的两个节点  
        if (pos1 >= 0)//第一个权值最小的节点  
        {  
          if (tempHuffmanTree[pos1].second < tempHuffmanTree[pos2].second)  
          {  
            min1i = pos1;  
            pos1--;  
          }  
          else  
          {  
            min1i = pos2;  
            pos2++;  
          }  
        }  
        else  
        {  
          min1i = pos2;  
          pos2++;  
        }  
        if (pos1 >= 0)//第二个权值最小的节点  
        {  

          if (tempHuffmanTree[pos1].second < tempHuffmanTree[pos2].second)  
          {  
            min2i = pos1;  
            pos1--;  
          }  
          else  
          {  
           min2i = pos2;  
           pos2++;  
         }  
       }  
       else  
       {  
        min2i = pos2;  
        pos2++; 
       }  
       tempHuffmanTree[len + i] = std::make_pair(std::to_string(len + i), tempHuffmanTree[min1].second + tempHuffmanTree[min2].second); 
       parent_node[min1i] = vocab_size + a;  
       parent_node[min2i] = vocab_size + a;  
       binary[min2i] = 1;//节点编码为1，之前默认是0。
       //add node to vocablist
       /*vocabWordObject tempVocabWordObject;
       tempVocabWordObject.word = std::to_string(len + i);
       tempVocabWordObject.num_occur = tempHuffmanTree[len + i];
       worker.add_object(vocabList, tempVocabWordObject);*/
    }

    // Now assign binary code to each vocabulary word  
    for (a = 0; a < len; a++)  
    {  
      b = a;  
      i = 0;  
      vocabWordObject* findObj = vocabList.find(tempHuffmanTree[a].first);
      while (1)  
      {  
        findObj.code[i] = binary[b];  
        findObj.point[i] = b;  
        i++;  
        b = parent_node[b];  
        if (b == len * 2 - 2) break;  
      }  

      findObj.codelen = i;  
      findObj.point[0] = len - 2;  
      for (b = 0; b < i; b++)  
      {  
        findObj.code[i - b - 1] = findObj.code[b];  
        findObj.point[i - b] = findObj.point[b] - len;  
      }  
    }  
}  

//这个线程函数执行之前，已经做好了一些工作：根据词频排序的词汇表，每个单词的huffman编码  
void TrainModelThread()  
{ 
  //========================================================主程序=======================================================//
  //step 01.创造exp_table的Aggregator
  std::vector<double> initExpTable;
  initExpTable.resize(EXP_TABLE_SIZE);

  for (i = 0; i <= EXP_TABLE_SIZE; i++)  
  {  
    //expTable[i] = exp((i -500)/ 500 * 6) 即 e^-6 ~ e^6  
    expTable[i] = exp((i / (double)EXP_TABLE_SIZE * 2 - 1) * MAX_EXP); // Precompute the exp() table  
    //expTable[i] = 1/(1+e^6) ~ 1/(1+e^-6)即 0.01 ~ 1 的样子  
    expTable[i] = expTable[i] / (expTable[i] + 1);                   // Precompute f(x) = x / (x + 1)  
  } 

  Husky::Aggregator<std::vector<double>> expTable(initExpTable,
    [](std::vector<double> &a, const std::vector<double> &b){
      for (int i = a.size() - 1; i >= 0; i--) a[i] += b[i];
    });

  //step 02.创造Vocab的Aggeragator
  std::vector<vocabWordObject> initVocab;
  initVocab.resize(sum.get_value());

  Husky::Aggregator<std::vector<vocabWordObject>> vocab(initVocab,
    [](std::vector<vocabWordObject> &a, const std::vector<vocabWordObject> &b) {
      a.push_back(b[0]);
    });

  worker.list_execute(vocabList, [&](vocabWordObject &tmpObj){
    std::vector<vocabWordObject> tempVec;
    tempVec.push_back(tmpObj);
    vocab.update(tempVec);
  });

  //step 03.创造Hash_table的Aggregator

  //step 04.创建其他向量的Aggregator
  real *neu1 = (real *)calloc(layer1_size, sizeof(real));  
  real *neu1e = (real *)calloc(layer1_size, sizeof(real));  
  val syn0Global =
      Array.fill[Float](vocabSize * vectorSize)((initRandom.nextFloat() - 0.5f) / vectorSize)
    val syn1Global = new Array[Float](vocabSize * vectorSize)

  //step 04.读文件创造固定长度的sentence

  //step 05.设置random的种子
  srand((unsigned int)time(NULL));


  //step 06.train cbow using HS and NG
  if (cbow)  
  {  //train the cbow architecture  
      // in -> hidden  
    if (hs){
      int pos = 0;
      while (pos < sentence.length) {
        int word = sentence[pos];
        int b = random.nextInt[window];
        // Train Skip-gram
        int a = b;

        //迭代相加
        while (a < window * 2 + 1 - b) {
          if (a != window) {
            int c = pos - window + a;
            if (c >= 0 && c < sentence.length) {
                int lastWord = sentence[c];
                if (lastWord != -1){  
                    for (c = 0; c < vectorSize; c++)//layer1_size词向量的维度，默认值是100  
                        neu1[c] += syn0[c + lastWord * vectorSize];//传说中的向量和？
                }
            }
          }
        }

        int d = 0;
        while (d < vocab.get_value()[word].codeLen) {
          int f = 0;
          int l2 = vocab.get_value()[word].point[d] * vectorSize;//point应该记录的是huffman的路径。找到当前节点，并算出偏移  
           // Propagate hidden -> output  
            for (c = 0; c < vectorSize; c++) f += neu1[c] * syn1[c + l2];//计算内积  
            if (f <= -MAX_EXP) continue;//内积不在范围内直接丢弃  
            else if (f >= MAX_EXP) continue;  
            else f = expTable.get_value()[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))];//内积之后sigmoid函数  
            // 'g' is the gradient multiplied by the learning rate  
            double g = (1 - vocab.get_value()[word].code[d] - f) * alpha;//偏导数的一部分  
            //layer1_size是向量的维度  
            // Propagate errors output -> hidden 反向传播误差，从huffman树传到隐藏层。下面就是把当前内节点的误差传播给隐藏层，syn1[c + l2]是偏导数的一部分。  
            for (c = 0; c < vectorSize; c++) neu1e[c] += g * syn1[c + l2];  
            // Learn weights hidden -> output 更新当前内节点的向量，后面的neu1[c]其实是偏导数的一部分  
            for (c = 0; c < vectorSize; c++) syn1[c + l2] += g * neu1[c];  
        }
    }

    // NEGATIVE SAMPLING  
    if (negative > 0){
      for (d = 0; d < negative + 1; d++)  
      {  
        if (d == 0)  
        {  
          target = word;//目标单词  
          label = 1;//正样本  
        }  
        else  
        {  
          next_random = next_random * (unsigned long long)25214903917 + 11;  
          target = table[(next_random >> 16) % table_size];  
          if (target == 0) target = next_random % (vocab_size - 1) + 1;  
          if (target == word) continue;  
          label = 0;//负样本  
        }  
        l2 = target * layer1_size;  
        f = 0;  
        for (c = 0; c < layer1_size; c++)  
            f += neu1[c] * syn1neg[c + l2];//内积  
        if (f > MAX_EXP)  
            g = (label - 1) * alpha;  
        else if (f < -MAX_EXP)  
            g = (label - 0) * alpha;  
        else g = (label - expTable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;  
        for (c = 0; c < layer1_size; c++)  
            neu1e[c] += g * syn1neg[c + l2];//隐藏层的误差  
        for (c = 0; c < layer1_size; c++)  
            syn1neg[c + l2] += g * neu1[c];//更新负样本向量  
      }  
      // hidden -> in  
      for (a = b; a < window * 2 + 1 - b; a++)  
      if (a != window)//cbow模型 更新的不是中间词语的向量，而是周围几个词语的向量。  
      {  
        c = sentence_position - window + a;  
        if (c < 0) continue;  
        if (c >= sentence_length) continue;  
        last_word = sen[c];  
        if (last_word == -1) continue;  
        for (c = 0; c < layer1_size; c++)  
            syn0[c + last_word * layer1_size] += neu1e[c];//更新词向量  
      }
    }  
    else  
    {  //train skip-gram  
       for (a = b; a < window * 2 + 1 - b; a++)  
       if (a != window)//扫描周围几个词语  
       {  
        c = sentence_position - window + a;  
        if (c < 0) continue;  
        if (c >= sentence_length) continue;  
        last_word = sen[c];  
        if (last_word == -1) continue;  
        l1 = last_word * layer1_size;  
        for (c = 0; c < layer1_size; c++)  
            neu1e[c] = 0;  
        // HIERARCHICAL SOFTMAX  
        if (hs)  
        for (d = 0; d < vocab[word].codelen; d++)//遍历叶子节点  
        {  
          f = 0;  
          l2 = vocab[word].point[d] * layer1_size;//point记录的是huffman的路径  
          // Propagate hidden -> output 感觉源代码这个英语注释有点误导人，这里的隐藏层就是输入层，就是词向量。  
          for (c = 0; c < layer1_size; c++)  
              f += syn0[c + l1] * syn1[c + l2];//计算两个词向量的内积  
          if (f <= -MAX_EXP) continue;  
          else if (f >= MAX_EXP) continue;  
          else f = expTable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))];  
          // 'g' is the gradient multiplied by the learning rate  
          g = (1 - vocab[word].code[d] - f) * alpha;//偏导数的一部分  
          // Propagate errors output -> hidden  
          for (c = 0; c < layer1_size; c++)  
              neu1e[c] += g * syn1[c + l2];//隐藏层的误差  
          // Learn weights hidden -> output  
          for (c = 0; c < layer1_size; c++)  
              syn1[c + l2] += g * syn0[c + l1];//更新叶子节点向量  
        }  
        // NEGATIVE SAMPLING  
        if (negative > 0)//这个同cobow差不多  
        for (d = 0; d < negative + 1; d++)  
        {  
          if (d == 0)  
          {  
            target = word;  
            label = 1;  
          }  
          else  
          {  
            next_random = next_random * (unsigned long long)25214903917 + 11;  
            target = table[(next_random >> 16) % table_size];  
            if (target == 0) target = next_random % (vocab_size - 1) + 1;  
            if (target == word) continue;  
            label = 0;  
          }  
          l2 = target * layer1_size;  
          f = 0;  
          for (c = 0; c < layer1_size; c++)  
              f += syn0[c + l1] * syn1neg[c + l2];  
          if (f > MAX_EXP) g = (label - 1) * alpha;  
          else if (f < -MAX_EXP)  
              g = (label - 0) * alpha;  
          else g = (label - expTable[(int)((f + MAX_EXP) * (EXP_TABLE_SIZE / MAX_EXP / 2))]) * alpha;  
          for (c = 0; c < layer1_size; c++)  
              neu1e[c] += g * syn1neg[c + l2];  
          for (c = 0; c < layer1_size; c++)  
              syn1neg[c + l2] += g * syn0[c + l1];  
        }  
  
        // Learn weights input -> hidden  
        for (c = 0; c < layer1_size; c++)  
            syn0[c + l1] += neu1e[c];//更新周围几个词语的向量  
      }  
    }  
    sentence_position++;  
    if (sentence_position >= sentence_length)  
    {  
      sentence_length = 0;  
      continue;  
    }  
  }  
  fclose(fi);  
  free(neu1);  
  free(neu1e);  
  pthread_exit(NULL);  
}  

void word2vec(){

    printf("WORD VECTOR estimation toolkit v 0.1b\n\n");  
    printf("Options:\n");  
    printf("Parameters for training:\n"); 

    //输出文件：词向量或者词聚类  
    printf("\t-output <file>\n");  
    printf("\t\tUse <file> to save the resulting word vectors / word clusters\n");  
  
    //=========================================================参数设置==================================================================//
    //词向量的维度，默认值是100  
    int wordVecSize = 100;
   /* printf("\t-size <int>\n");  
    printf("\t\tSet size of word vectors; default is 100\n");  */
  
    //窗口大小，默认是5  
    int windowSize = 5;
    /*printf("\t-window <int>\n");  
    printf("\t\tSet max skip length between words; default is 5\n");  */
  
    //设定词出现频率的阈值，对于常出现的词会被随机下采样  
    double threshold = 1e-5;
    /*printf("\t-sample <float>\n");  
    printf("\t\tSet threshold for occurrence of words. Those that appear with higher frequency");  
    printf(" in the training data will be randomly down-sampled; default is 0 (off), useful value is 1e-5\n");  */
  
    //是否采用softmax体系  
    int isSoftmax = 0;
    /*printf("\t-hs <int>\n");  
    printf("\t\tUse Hierarchical Softmax; default is 1 (0 = not used)\n");  */
  
    //负样本的数量，默认是0，通常使用5-10。0表示不使用。 
    int negSampleNum = 0; 
    /*printf("\t-negative <int>\n");  
    printf("\t\tNumber of negative examples; default is 0, common values are 5 - 10 (0 = not used)\n"); */ 
  
    //开启的线程数量  
    /*printf("\t-threads <int>\n");  
    printf("\t\tUse <int> threads (default 1)\n");  */
    //最小阈值。对于出现次数少于该值的词，会被抛弃掉。  
    int minThreshold = 5;
   /* printf("\t-min-count <int>\n");  
    printf("\t\tThis will discard words that appear less than <int> times; default is 5\n");  */
  
    //学习速率初始值，默认是0.025  
    double alpha = 0.025;
    /*printf("\t-alpha <float>\n");  
    printf("\t\tSet the starting learning rate; default is 0.025\n");  */
  
    //输出词类别，而不是词向量  
    int wordCatagory = 0;
    /*printf("\t-classes <int>\n");  
    printf("\t\tOutput word classes rather than word vectors; default number of classes is 0 (vectors are written)\n");  */
  
    //debug模式，默认是2，表示在训练过程中会输出更多信息  
    int debug = 2;
    /*printf("\t-debug <int>\n");  
    printf("\t\tSet the debug mode (default = 2 = more info during training)\n");  */
  
    //是否用binary模式保存数据，默认是0，表示否。  
    int binary = 0;
    /*printf("\t-binary <int>\n");  
    printf("\t\tSave the resulting vectors in binary moded; default is 0 (off)\n");  */
  
    //保存词汇到这个文件  
    std::string save_vocab_file = "output.txt";
    /*printf("\t-save-vocab <file>\n");  
    printf("\t\tThe vocabulary will be saved to <file>\n");  */
  
    //词汇从该文件读取，而不是由训练数据重组  
    std::string read_vocab_file = "input.txt";
    /*printf("\t-read-vocab <file>\n");  
    printf("\t\tThe vocabulary will be read from <file>, not constructed from the training data\n");  */
  
    //是否采用continuous bag of words算法。默认是0，表示采用另一个叫skip-gram的算法。  
    int isCBOM = 0;
    int isSkipGram = 1;
    /*printf("\t-cbow <int>\n");  
    printf("\t\tUse the continuous bag of words model; default is 0 (skip-gram model)\n");  */
  
    //工具使用样例  
    /*printf("\nExamples:\n");  
    printf("./word2vec -train data.txt -output vec.txt -debug 2 -size 200 -window 5 -sample 1e-4 -negative 5 -hs 0 -binary 0 -cbow 1\n\n");  
    return 0;  */

    //=====================================================变量定义========================================================//
    /*char train_file[MAX_STRING], output_file[MAX_STRING];  
    char save_vocab_file[MAX_STRING], read_vocab_file[MAX_STRING];  
    struct vocab_word *vocab;  //对于语料库中的每一个词放一个词的struct
    int binary = 0, cbow = 0, debug_mode = 2, window = 5, min_count = 5, num_threads = 1, min_reduce = 1;  
    int *vocab_hash;  
    long long vocab_max_size = 1000, vocab_size = 0, layer1_size = 100;  
    long long train_words = 0, word_count_actual = 0, file_size = 0, classes = 0;  
    real alpha = 0.025, starting_alpha, sample = 0;  
    real *syn0, *syn1, *syn1neg, *expTable;  
    clock_t start;  

    int hs = 1, negative = 0; 
    const int table_size = 1e8;  
    int *table;  */

    //==================================================链表初始化=================================================//
    /*vocab = (struct vocab_word *)calloc(vocab_max_size, sizeof(struct vocab_word));  
    vocab_hash = (int *)calloc(vocab_hash_size, sizeof(int));  
    expTable = (real *)malloc((EXP_TABLE_SIZE + 1) * sizeof(real));  
    for (i = 0; i < EXP_TABLE_SIZE; i++)  
    {  
    //expTable[i] = exp((i -500)/ 500 * 6) 即 e^-6 ~ e^6  
    expTable[i] = exp((i / (real)EXP_TABLE_SIZE * 2 - 1) * MAX_EXP); // Precompute the exp() table  
    //expTable[i] = 1/(1+e^6) ~ 1/(1+e^-6)即 0.01 ~ 1 的样子  
    expTable[i] = expTable[i] / (expTable[i] + 1);                   // Precompute f(x) = x / (x + 1)  
    }  *//***************************斟酌一下这里咋写*/
  //TrainModel();  
    //=================================================主程序======================================================//

    //step 1.从训练文件学习词汇,主要制作词汇表学习词频

    Husky::BaseWorker & worker = Husky::Context::get_worker<Husky::BaseWorker>();
    auto & vocabList = worker.create_list<vocabWordObject>("vocab_list");
    LearnVocabFromTrainFile(vocabList);

  InitNet();  
  if (negative > 0) InitUnigramTable();  
  start = clock();  
  for (a = 0; a < num_threads; a++) pthread_create(&pt[a], NULL, TrainModelThread, (void *)a);  
  for (a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);  
  fo = fopen(output_file, "wb");  
  if (classes == 0) //不需要聚类，只需要输出词向量  
  {  
    // Save the word vectors  
    fprintf(fo, "%lld %lld\n", vocab_size, layer1_size);  
    for (a = 0; a < vocab_size; a++)  
    {  
      fprintf(fo, "%s ", vocab[a].word);  
      if (binary)  
          for (b = 0; b < layer1_size; b++)  
              fwrite(&syn0[a * layer1_size + b], sizeof(real), 1, fo);  
      else  
          for (b = 0; b < layer1_size; b++)  
              fprintf(fo, "%lf ", syn0[a * layer1_size + b]);  
      fprintf(fo, "\n");  
    }  
  }  
  else //使用k-means进行聚类  
  {  
    // Run K-means on the word vectors  
    int clcn = classes, iter = 10, closeid;  
    int *centcn = (int *)malloc(classes * sizeof(int));//该类别的数量  
    int *cl = (int *)calloc(vocab_size, sizeof(int));//词到类别的映射  
    real closev, x;  
    real *cent = (real *)calloc(classes * layer1_size, sizeof(real));//质心数组  
    for (a = 0; a < vocab_size; a++)  
        cl[a] = a % clcn;//任意分类？  
    for (a = 0; a < iter; a++)  
    {  
      for (b = 0; b < clcn * layer1_size; b++)  
          cent[b] = 0;//质心清零  
      for (b = 0; b < clcn; b++)  
          centcn[b] = 1;  
      for (c = 0; c < vocab_size; c++)  
      {  
        for (d = 0; d < layer1_size; d++)  
            cent[layer1_size * cl[c] + d] += syn0[c * layer1_size + d];//求和放到质心数组中  
        centcn[cl[c]]++;//类别数量加1  
      }  
      for (b = 0; b < clcn; b++)//遍历所有类别  
      {  
        closev = 0;  
        for (c = 0; c < layer1_size; c++)  
        {  
          cent[layer1_size * b + c] /= centcn[b];//均值，就是求新的质心  
          closev += cent[layer1_size * b + c] * cent[layer1_size * b + c];  
        }  
        closev = sqrt(closev);  
        for (c = 0; c < layer1_size; c++)  
            cent[layer1_size * b + c] /= closev;//对质心进行归一化？  
      }  
      for (c = 0; c < vocab_size; c++)//对所有词语重新分类  
      {  
        closev = -10;  
        closeid = 0;  
        for (d = 0; d < clcn; d++)  
        {  
          x = 0;  
          for (b = 0; b < layer1_size; b++)  
              x += cent[layer1_size * d + b] * syn0[c * layer1_size + b];//内积  
          if (x > closev)  
          {  
            closev = x;  
            closeid = d;  
          }  
        }  
        cl[c] = closeid;  
      }  
    }  
    // Save the K-means classes  
    for (a = 0; a < vocab_size; a++)  
        fprintf(fo, "%s %d\n", vocab[a].word, cl[a]);  
    free(centcn);  
    free(cent);  
    free(cl);  
  }  
  fclose(fo);  
}  

int main(int argc, char ** argv) {
    Husky::run_job(word2vec, argv[1]);
    return 0;
}