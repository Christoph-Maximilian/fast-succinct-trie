#include <FST.hpp>
#include <bitset>
#include <fstream>
#include <sstream>

FST::FST(int cutoff_level) : cutoff_level_(cutoff_level), nodeCountU_(0), childCountU_(0),
                             cbitsU_(NULL), tbitsU_(NULL), obitsU_(NULL), valuesU_(NULL),
                             cbytes_(NULL), tbits_(NULL), sbits_(NULL), values_(NULL),
                             tree_height_(0), last_value_pos_(0),
                             c_lenU_(0), o_lenU_(0), c_memU_(0), t_memU_(0), o_memU_(0), val_memU_(0),
                             c_mem_(0), t_mem_(0), s_mem_(0), val_mem_(0), num_t_(0), number_values(0) {}

FST::~FST() {
    if (cbitsU_) delete cbitsU_;
    if (tbitsU_) delete tbitsU_;
    if (obitsU_) delete obitsU_;
    if (valuesU_) delete valuesU_;

    if (cbytes_) delete cbytes_;
    if (tbits_) delete tbits_;
    if (sbits_) delete sbits_;
    if (values_) delete values_;
}

//stat
uint64_t FST::cMemU() { return c_memU_; }

uint64_t FST::tMemU() { return t_memU_; }

uint64_t FST::oMemU() { return o_memU_; }

uint64_t FST::keyMemU() { return (c_memU_ + t_memU_ + o_memU_); }

uint64_t FST::valueMemU() { return static_cast<uint32_t >(values_U_succinct.size() * 8); }

uint64_t FST::cMem() { return c_mem_; }

uint64_t FST::tMem() { return t_mem_; }

uint64_t FST::sMem() { return s_mem_; }

uint64_t FST::keyMem() { return (c_mem_ + t_mem_ + s_mem_); }

uint64_t FST::valueMem() { return val_mem_; }

uint64_t FST::mem() { return (c_memU_ + t_memU_ + o_memU_ + val_memU_ + c_mem_ + t_mem_ + s_mem_ + val_mem_); }

uint32_t FST::numT() { return num_t_; }


//*******************************************************************
inline bool
FST::insertChar_cond(const uint8_t ch, vector<uint8_t> &c, vector<uint64_t> &t, vector<uint64_t> &s, uint64_t &pos,
                     uint64_t &nc) {
    if (c.empty() || c.back() != ch) {
        c.push_back(ch);
        if (c.size() == 1) {
            setBit(s.back(), pos % 64);
            nc++;
        }
        pos++;
        if (pos % 64 == 0) {
            t.push_back(0);
            s.push_back(0);
        }
        return true;
    } else {
        if (pos % 64 == 0)
            setBit(t.rbegin()[1], 63);
        else
            setBit(t.back(), (pos - 1) % 64);
        return false;
    }
}

inline bool
FST::insertChar(const uint8_t ch, bool isTerm, vector<uint8_t> &c, vector<uint64_t> &t, vector<uint64_t> &s,
                uint64_t &pos, uint64_t &nc, bool set_SBit) {
    c.push_back(ch);
    if (!isTerm)
        setBit(t.back(), pos % 64);
    if (set_SBit) {
        setBit(s.back(), pos % 64);
        nc++;
    }

    pos++;
    if (pos % 64 == 0) {
        t.push_back(0);
        s.push_back(0);
    }
    return true;
}

//******************************************************
// STATISTICS
//******************************************************

std::string FST::export_stats() {
    std::stringstream stream;
    stream << "mem = " << this->mem() << std::endl;

    // DENSE NODES
    stream << "cMemU = " << this->cMemU() << std::endl;
    stream << "tMemU = " << this->tMemU() << std::endl;
    stream << "oMemU = " << this->oMemU() << std::endl;
    stream << "keyMemU = " << this->keyMemU() << std::endl;
    stream << "valueMemU = " << this->valueMemU() << std::endl;

    // SPARSE NODES
    stream << "cMem = " << this->cMem() << std::endl;
    stream << "tMem = " << this->tMem() << std::endl;
    stream << "sMem = " << this->sMem() << std::endl;
    stream << "keyMem = " << this->keyMem() << std::endl;
    stream << "valueMem = " << this->valueMem() << std::endl;
    return stream.str();
}

//******************************************************
// LOAD
//******************************************************

uint8_t denorm_levels[4] = {1, 4, 16, 64};

void FST::load(vector<uint8_t> &keys, vector<uint64_t> &values, int longestKeyLen) {
    number_values = values.size();
    tree_height_ = static_cast<uint32_t >(longestKeyLen);
    vector<vector<uint8_t> > c;
    vector<vector<uint64_t> > t;
    vector<vector<uint64_t> > s;
    vector<vector<uint64_t> > val;

    vector<std::pair<vector<uint64_t>, vector<uint64_t >>> val_succ_tmp(longestKeyLen);
    vector<uint64_t> values_count_per_level(longestKeyLen);

    vector<uint64_t> keys_per_level;
    vector<uint64_t> pos_list;
    vector<uint64_t> nc; //node count

    // init
    for (int i = 0; i < longestKeyLen; i++) {
        c.emplace_back(vector<uint8_t>());
        t.emplace_back(vector<uint64_t>());
        s.emplace_back(vector<uint64_t>());
        val.emplace_back(vector<uint64_t>());
        keys_per_level.push_back(0);

        pos_list.push_back(0);
        nc.push_back(0);

        t[i].push_back(0);
        s[i].push_back(0);
    }

    int last_value_level = 0;
    uint64_t index = 0;

    uint8_t key[8];
    uint64_t value_index = std::numeric_limits<uint64_t>::max();
    while (index < keys.size()) {
        value_index++;
        auto info_byte = keys[index];
        auto number_denorm_cells = denorm_levels[info_byte >> 6];
        auto number_common_cells = info_byte & 15;
        auto key_size = number_common_cells + 1;

        for (auto ncl = 0; ncl < number_common_cells; ncl++) {
            index++;
            key[ncl] = keys[index];
        }

        index++;
        uint64_t end_index = index + static_cast<uint64_t>(number_denorm_cells);


        // set denormalize cell at the end
        key[number_common_cells] = keys[index];
        uint64_t value = values[value_index];

        int i = 0;
        while (i < key_size && !insertChar_cond((uint8_t) key[i], c[i], t[i], s[i], pos_list[i], nc[i])) {
            i++;
        }

        if (i < key_size) {
            //  todo commented out: if (k + 1 < (int) key_size) {...
            auto cpl = static_cast<uint32_t >(key_size); //commonPrefixLen(key, keys[k + 1]);
            if (i < (cpl - 1)) {
                // There will be a child node -> set T bit of corresponding position within node
                if (pos_list[i] % 64 == 0)
                    setBit(t[i].rbegin()[1], 63);
                else
                    setBit(t[i].back(), (pos_list[i] - 1) % 64);
            }
            // the first (denormalized) key MAY be inserted in insertChar
            if (i == key_size - 1) {
                index++;
            }
            while (i < (cpl - 1)) {
                i++;
                if (i < (cpl - 1))
                    insertChar((uint8_t) key[i], false, c[i], t[i], s[i], pos_list[i], nc[i], true);
            }
            for (; index < end_index; index++) {
                key[number_common_cells] = keys[index];
                insertChar((uint8_t) key[i], true, c[i], t[i], s[i], pos_list[i], nc[i],
                           end_index - index == number_denorm_cells);
            }

            // todo: save values more efficiently here
            for (auto nv = 0; nv < 1; nv++) {
                values_count_per_level[i] += number_denorm_cells;
                if (val_succ_tmp[i].first.empty()) {
                    val_succ_tmp[i].first.push_back(value);
                    val_succ_tmp[i].second.push_back(number_denorm_cells);
                    break;
                }
                if (*val_succ_tmp[i].first.rbegin() == value) {
                    (*val_succ_tmp[i].second.rbegin()) += number_denorm_cells;
                } else {
                    val_succ_tmp[i].first.push_back(value);
                    val_succ_tmp[i].second.push_back(number_denorm_cells);
                }
            }
        } else
            cout << "ERROR!\n";

        if (index >= keys.size() - 1)
            last_value_level = i;

    }

    // put together
    uint64_t nc_total = 0;
    for (int i = 0; i < (int) nc.size(); i++)
        nc_total += nc[i];

    cout << "cutoff_level_ = " << cutoff_level_ << "\n";

    // determine the position of the last value for range query boundary check
    if (last_value_level < cutoff_level_) {
        for (int i = 0; i <= last_value_level; i++)
            last_value_pos_ -= val[i].size();
        last_value_pos_++;
    } else {
        for (int i = cutoff_level_; i <= last_value_level; i++)
            last_value_pos_ += val[i].size();
        last_value_pos_--;
    }

    //-------------------------------------------------
    vector<vector<uint64_t> > cU;
    vector<vector<uint64_t> > tU;
    vector<vector<uint8_t> > oU;

    for (int i = 0; i < cutoff_level_; i++) {
        cU.push_back(vector<uint64_t>());
        tU.push_back(vector<uint64_t>());
        oU.push_back(vector<uint8_t>());
    }

    for (int i = 0; i < cutoff_level_; i++) {
        uint64_t sbitPos_i = 0;
        uint64_t smaxPos_i = pos_list[i];

        uint8_t ch = 0;
        uint64_t j = 0;
        uint64_t k = 0;
        for (j = 0; j < s[i].size(); j++) {
            if (sbitPos_i >= smaxPos_i) break;
            for (k = 0; k < 64; k++) {
                if (sbitPos_i >= smaxPos_i) break;

                ch = c[i][sbitPos_i];
                if (readBit(s[i][j], k)) {
                    for (int l = 0; l < 4; l++) {
                        cU[i].push_back(0);
                        tU[i].push_back(0);
                    }
                    nodeCountU_++;
                    oU[i].push_back(0);
                    setLabel(cU[i].data() + cU[i].size() - 4, ch);
                    if (readBit(t[i][j], k)) {
                        setLabel(tU[i].data() + tU[i].size() - 4, ch);
                        childCountU_++;
                    }
                } else {
                    setLabel(cU[i].data() + cU[i].size() - 4, ch);
                    if (readBit(t[i][j], k)) {
                        setLabel(tU[i].data() + tU[i].size() - 4, ch);
                        childCountU_++;
                    }
                }
                sbitPos_i++;
            }
        }
    }

    o_lenU_ = 0;
    u_int64_t vallenU = 0;
    for (uint64_t i = 0; i < cU.size(); i++) {
        c_lenU_ += cU[i].size();
        o_lenU_ += oU[i].size();
        vallenU += values_count_per_level[i];
    }

    auto *cbitsU = new uint64_t[c_lenU_];
    auto *tbitsU = new uint64_t[c_lenU_];
    auto *obitsU = new uint64_t[o_lenU_ / 64 + 1];
    vbitsU_bits.resize(vallenU / 64 + 1, 0);

    // init
    for (uint64_t i = 0; i < c_lenU_; i++) {
        cbitsU[i] = 0;
        tbitsU[i] = 0;
    }
    for (uint64_t i = 0; i < vallenU / 64 + 1; i++) {
        vbitsU_bits[i] = 0;
    }


    for (uint64_t i = 0; i < (o_lenU_ / 64 + 1); i++)
        obitsU[i] = 0;

    uint64_t c_bitPosU = 0;
    for (uint64_t i = 0; i < cU.size(); i++) {
        for (uint64_t j = 0; j < cU[i].size(); j++) {
            for (uint64_t k = 0; k < 64; k++) {
                if (readBit(cU[i][j], k)) {
                    setBit(cbitsU[c_bitPosU / 64], c_bitPosU % 64);
                }
                c_bitPosU++;
            }
        }
    }

    u_int64_t c_sizeU = (c_lenU_ / 32 + 1) * 32; // round-up to 1024-bit block size for Poppy
    cbitsU_ = new BitmapRankFPoppy(cbitsU, c_sizeU * 64);
    c_memU_ = cbitsU_->getNbits() / 8; //stat

    uint64_t t_bitPosU = 0;
    for (uint64_t i = 0; i < tU.size(); i++) {
        for (uint64_t j = 0; j < tU[i].size(); j++) {
            for (uint64_t k = 0; k < 64; k++) {
                if (readBit(tU[i][j], k))
                    setBit(tbitsU[t_bitPosU / 64], t_bitPosU % 64);
                t_bitPosU++;
            }
        }
    }

    assert (t_bitPosU < c_bitPosU);
    u_int64_t t_sizeU = (c_lenU_ / 32 + 1) * 32; // round-up to 1024-bit block size for Poppy
    tbitsU_ = new BitmapRankFPoppy(tbitsU, t_sizeU * 64);
    t_memU_ = tbitsU_->getNbits() / 8; //stat

    //todo remove this part?
    uint64_t o_bitPosU = 0;
    for (uint64_t i = 0; i < oU.size(); i++) {
        for (uint64_t j = 0; j < oU[i].size(); j++) {
            if (oU[i][j] == 1)
                setBit(obitsU[o_bitPosU / 64], o_bitPosU % 64);
            o_bitPosU++;
        }
    }

    u_int64_t o_sizeU = (o_lenU_ / 64 / 32 + 1) * 32; // round-up to 1024-bit block size for Poppy
    obitsU_ = new BitmapRankFPoppy(obitsU, o_sizeU * 64);
    o_memU_ = obitsU_->getNbits() / 8; //stat

    uint64_t val_posU = 0;
    uint64_t bit_value_posU = 0;
    uint64_t last_value = 0;


    for (u_int32_t i = 0; i < cutoff_level_; i++) {
        //concat values
        values_U_succinct.insert(values_U_succinct.end(), val_succ_tmp[i].first.begin(), val_succ_tmp[i].first.end());
        // set bitmap
        for (uint64_t j = 0; j < val_succ_tmp[i].second.size(); j++) {
            bit_value_posU += val_succ_tmp[i].second[j];
            setBit(vbitsU_bits[bit_value_posU / 64], bit_value_posU % 64);
        }
    }

    u_int64_t v_sizeU = (bit_value_posU / 64 / 32 + 1) * 32;
    vbitsU_bits.resize(v_sizeU, 0);
    vbitsU_ = new BitmapRankFPoppy(vbitsU_bits.data(), v_sizeU * 64);

    val_memU_ = static_cast<uint32_t >(values_U_succinct.size() * 8 + vbitsU_->getNbits() / 8);

    //-------------------------------------------------
    for (int i = cutoff_level_; i < c.size(); i++)
        c_mem_ += c[i].size();

    for (int i = cutoff_level_; i < c.size(); i++)
        val_mem_ += values_count_per_level[i] * sizeof(uint64_t);

    if (c_mem_ % 64 == 0) {
        t_mem_ = c_mem_ / 64;
        s_mem_ = c_mem_ / 64;
    } else {
        t_mem_ = c_mem_ / 64 + 1;
        s_mem_ = c_mem_ / 64 + 1;
    }

    t_mem_ = (t_mem_ / 32 + 1) * 32; // round-up to 2048-bit block size for Poppy
    s_mem_ = (s_mem_ / 32 + 1) * 32; // round-up to 2048-bit block size for Poppy

    cbytes_ = new uint8_t[c_mem_];
    auto *tbits = new uint64_t[t_mem_];
    auto *sbits = new uint64_t[s_mem_];
    values_ = new uint64_t[val_mem_ / sizeof(uint64_t)];

    // init
    for (uint64_t i = 0; i < t_mem_; i++) {
        tbits[i] = 0;
        sbits[i] = 0;
    }

    //todo: parallelize the following stuff

    uint64_t c_pos = 0;
    for (int i = cutoff_level_; i < (int) c.size(); i++) {
        for (uint64_t j = 0; j < c[i].size(); j++) {
            cbytes_[c_pos] = c[i][j];
            c_pos++;
        }
    }

    uint64_t t_bitPos = 0;
    for (int i = cutoff_level_; i < (int) t.size(); i++) {
        uint64_t bitPos_i = 0;
        uint64_t maxPos_i = pos_list[i];
        for (uint64_t j = 0; j < t[i].size(); j++) {
            if (bitPos_i >= maxPos_i) break;
            for (uint64_t k = 0; k < 64; k++) {
                if (bitPos_i >= maxPos_i) break;
                if (readBit(t[i][j], k))
                    setBit(tbits[t_bitPos / 64], t_bitPos % 64);
                t_bitPos++;
                bitPos_i++;
            }
        }
    }

    tbits_ = new BitmapRankPoppy(tbits, t_mem_ * 64);
    t_mem_ = tbits_->getMem(); //stat

    uint64_t s_bitPos = 0;
    for (int i = cutoff_level_; i < (int) s.size(); i++) {
        uint64_t bitPos_i = 0;
        uint64_t maxPos_i = pos_list[i];
        for (uint64_t j = 0; j < s[i].size(); j++) {
            if (bitPos_i >= maxPos_i) break;
            for (int k = 0; k < 64; k++) {
                if (bitPos_i >= maxPos_i) break;
                if (readBit(s[i][j], k))
                    setBit(sbits[s_bitPos / 64], s_bitPos % 64);
                s_bitPos++;
                bitPos_i++;
            }
        }
    }

    sbits_ = new BitmapSelectPoppy(sbits, s_mem_ * 64);
    s_mem_ = sbits_->getMem(); //stat

    uint64_t val_pos = 0;
    for (int i = cutoff_level_; i < (int) val.size(); i++) {
        for (int j = 0; j < (int) val_succ_tmp[i].second.size(); j++) {
            for (uint64_t k = 0; k < val_succ_tmp[i].second[j]; k++) {
                values_[val_pos] = val_succ_tmp[i].first[j];
                val_pos++;
            }
        }
    }
    //-------------------------------------------------
    //node counts per level
    nc_ = nc;
    for (auto i = 0; i < c.size(); i++) {
        keys_per_level[i] = c[i].size();
    }
    keys_per_level_ = keys_per_level;

}

void FST::load(vector<uint64_t> &keys, vector<uint64_t> &values) {
    vector<string> keys_str;
    for (int i = 0; i < (int) keys.size(); i++) {
        char key[8];
        reinterpret_cast<uint64_t *>(key)[0] = __builtin_bswap64(keys[i]);
        keys_str.push_back(string(key, 8));
    }
    //load(keys_str, values, sizeof(uint64_t));
}

//******************************************************
// IS C BIT SET U(pper Level)? -> dense
//******************************************************
inline bool FST::isCbitSetU(uint64_t nodeNum, uint8_t kc) {
    return isLabelExist(cbitsU_->bits_ + (nodeNum << 2), kc);
}

//******************************************************
// IS T BIT SET U? dense
//******************************************************
inline bool FST::isTbitSetU(uint64_t nodeNum, uint8_t kc) {
    return isLabelExist(tbitsU_->bits_ + (nodeNum << 2), kc);
}

//******************************************************
// IS O BIT SET U? dense
//******************************************************
inline bool FST::isObitSetU(uint64_t nodeNum) {
    return readBit(obitsU_->bits_[nodeNum >> 6], nodeNum & (uint64_t) 63);
}

//******************************************************
// IS S BIT SET? sparse
//******************************************************
inline bool FST::isSbitSet(uint64_t pos) {
    return readBit(sbits_->bits_[pos >> 6], pos & (uint64_t) 63);
}

//******************************************************
// IS T BIT SET? sparse
//******************************************************
inline bool FST::isTbitSet(uint64_t pos) {
    return readBit(tbits_->bits_[pos >> 6], pos & (uint64_t) 63);
}

//******************************************************
// GET VALUE POS U dense
//******************************************************
inline uint64_t FST::valuePosU(uint64_t nodeNum, uint64_t pos) {
    return vbitsU_->rank(cbitsU_->rank(pos + 1)-tbitsU_->rank(pos + 1));
}

//******************************************************
// GET VALUE POS sparse
//******************************************************
inline uint64_t FST::valuePos(uint64_t pos) {
    return pos - tbits_->rank(pos + 1);
}

//******************************************************
// CHILD NODE NUM dense
//******************************************************
inline uint64_t FST::childNodeNumU(uint64_t pos) {
    return tbitsU_->rank(pos + 1);
}

inline uint64_t FST::childNodeNum(uint64_t pos) {
    return tbits_->rank(pos + 1);
}

//******************************************************
// CHILD POS sparse
//******************************************************
inline uint64_t FST::childpos(uint64_t nodeNum) {
    if (nodeNum - nodeCountU_ + 1 < 0) {
        std::cout << "select of negative value: nodenum=" << nodeNum << "nodecountU=" << nodeCountU_ << std::endl;
    }
    return sbits_->select(nodeNum - nodeCountU_ + 1);
}


//******************************************************
// NODE SIZE
//******************************************************
inline int FST::nodeSize(uint64_t pos) {
    pos++;
    //shift 6 to right, because 2^6 = 64 == 1 uint64_t,
    uint64_t startIdx = pos >> 6;
    // for the shift, only the six lowest level bits are relevant
    uint64_t shift = pos & (uint64_t) 0x3F;
    if (startIdx > sbits_->nbits_ / sbits_->kWordSize) {
        std::cout << "ERROR sbits pos=" << pos << " startIdx=" << startIdx << "size of s lut: " <<  sbits_->nbits_ / sbits_->kWordSize + 1 << std::endl;
    }
    uint64_t bits = sbits_->bits_[startIdx] << shift;

    if (bits > 0) //counting the number of leading zeros + 1 = nodesize
        return __builtin_clzll(bits) + 1;

    for (int i = 1; i < 5; i++) {
        bits = sbits_->bits_[startIdx + i];
        if (bits > 0)
            return 64 * i - shift + __builtin_clzll(bits) + 1;
    }
    return -1;
}

//******************************************************
// SIMD SEARCH
//******************************************************
inline bool FST::simdSearch(uint64_t &pos, uint64_t size, uint8_t target) {
    uint64_t s = 0;
    while (size >> 4) {
        __m128i cmp = _mm_cmpeq_epi8(_mm_set1_epi8(target),
                                     _mm_loadu_si128(reinterpret_cast<__m128i *>(cbytes_ + pos + s)));
        unsigned bitfield = _mm_movemask_epi8(cmp);
        if (bitfield) {
            pos += (s + __builtin_ctz(bitfield));
            return true;
        }
        s += 16;
        size -= 16;
    }

    if (size > 0) {
        __m128i cmp = _mm_cmpeq_epi8(_mm_set1_epi8(target),
                                     _mm_loadu_si128(reinterpret_cast<__m128i *>(cbytes_ + pos + s)));
        unsigned bitfield = _mm_movemask_epi8(cmp) & ((1 << size) - 1);
        if (bitfield) {
            pos += (s + __builtin_ctz(bitfield));
            return true;
        }
    }
    return false;
}

//******************************************************
// BINARY SEARCH
//******************************************************
inline bool FST::binarySearch(uint64_t &pos, uint64_t size, uint8_t target) {
    int64_t l = pos;
    int64_t r = pos + size - 1;
    int64_t m = (l + r) >> 1;

    while (l <= r) {
        if (cbytes_[m] == target) {
            pos = m;
            return true;
        } else if (cbytes_[m] < target)
            l = m + 1;
        else
            r = m - 1;
        m = (l + r) >> 1;
    }
    return false;
}

inline bool FST::binarySearch_lowerBound(uint64_t &pos, uint64_t size, uint8_t target) {
    int64_t rightBound = pos + size;
    int64_t l = pos;
    int64_t r = pos + size - 1;
    int64_t m = (l + r) >> 1;

    while (l < r) {
        if (cbytes_[m] == target) {
            pos = m;
            return true;
        } else if (cbytes_[m] < target)
            l = m + 1;
        else
            r = m - 1;
        m = (l + r) >> 1;
    }

    if (cbytes_[m] < target)
        pos = m + 1;
    else
        pos = m;

    return pos < rightBound;
}

//******************************************************
// LINEAR SEARCH
//******************************************************
inline bool FST::linearSearch(uint64_t &pos, uint64_t size, uint8_t target) {
    for (int i = 0; i < size; i++) {
        if (cbytes_[pos] == target)
            return true;
        pos++;
    }
    return false;
}

inline bool FST::linearSearch_lowerBound(uint64_t &pos, uint64_t size, uint8_t target) {
    for (int i = 0; i < size; i++) {
        if (cbytes_[pos] >= target)
            return true;
        pos++;
    }
    return false;
}

//******************************************************
// NODE SEARCH
//******************************************************
inline bool FST::nodeSearch(uint64_t &pos, int size, uint8_t target) {
    //Todo: Fix the binary search algorithm
    if (size < 7)
        return linearSearch(pos, size, target);
    else if (size < 12)
        return binarySearch(pos, size, target);
    else
        return simdSearch(pos, size, target);
}

inline bool FST::nodeSearch_lowerBound(uint64_t &pos, int size, uint8_t target) {
    if (size < 3)
        return linearSearch_lowerBound(pos, size, target);
    else
        return binarySearch_lowerBound(pos, size, target);
}


//******************************************************
// LOOKUP
//******************************************************
bool FST::lookup(const uint8_t *key, const int keylen, uint64_t &value) {
    int keypos = 0;

    uint64_t nodeNum = 0;
    auto kc = (uint8_t) key[keypos];
    uint64_t pos = kc;

    //******************************************************
    // SEARCH IN DENSE NODES
    // ******************************************************
    while (keypos < keylen && keypos < cutoff_level_) {
        kc = (uint8_t) key[keypos];
        pos = (nodeNum << 8) + kc;

        __builtin_prefetch(tbitsU_->bits_ + (nodeNum << 2) + (kc >> 6), 0);
        __builtin_prefetch(tbitsU_->rankLUT_ + ((pos + 1) >> 6), 0);

        if (!isCbitSetU(nodeNum, kc)) { // is C-label set?
            return false;
        }
        if (!isTbitSetU(nodeNum, kc)) {
            // this returns true if there is no child -> value can be found
            value = values_U_succinct[valuePosU(nodeNum, pos)];
            return true;
        }

        nodeNum = childNodeNumU(pos);
        keypos++;
    }

    //******************************************************
    // SEARCH IN SPARSE NODES
    // ******************************************************

    pos = (cutoff_level_ == 0) ? 0 : childpos(nodeNum);
    uint64_t tmp_pos = pos;
    while (keypos < keylen) {

        if (pos > 100000000000) {
            std::cout << "ERROR current pos=" << pos << " previous pos=" << tmp_pos << " childcountU=" << childCountU_
            << "keyposition=" << keypos << std::endl;
        }
        kc = (uint8_t) key[keypos];
        int nsize = nodeSize(pos);

        if (!nodeSearch(pos, nsize, kc)) {
            return false;
        }

        if (!isTbitSet(pos)) {
            //todo careful - we read out of bitvector here
            //uint64_t value_index = valuePos(pos);
            //assert(value_index < number_values);
            value = values_[valuePos(pos)];

            return true;
        }

        tmp_pos = pos;
        pos = childpos(childNodeNum(pos) + childCountU_);
        keypos++;

        __builtin_prefetch(cbytes_ + pos, 0, 1);
        __builtin_prefetch(tbits_->bits_ + (pos >> 6), 0, 1);
        __builtin_prefetch(tbits_->rankLUT_ + ((pos + 1) >> 9), 0);
    }

    //TODO: this is the important part -> here we detect that , added by @Christoph
    if (!isTbitSet(pos)) {
        //uint64_t value_index = valuePos(pos);
        //assert(value_index < number_values);
        value = values_[valuePos(pos)];
        return true;
    }
    return false;
}

bool FST::lookup(const uint64_t key, uint64_t &value) {
    uint8_t key_str[8];
    reinterpret_cast<uint64_t *>(key_str)[0] = __builtin_bswap64(key);

    return lookup(key_str, 8, value);
}


//******************************************************
// NEXT ITEM U
//******************************************************
inline bool FST::nextItemU(uint64_t nodeNum, uint8_t kc, uint8_t &cc) {
    return isLabelExist_lowerBound(cbitsU_->bits_ + (nodeNum << 2), kc, cc);
}

//******************************************************
// NEXT LEFT ITEM
//******************************************************
inline bool FST::nextLeftU(int level, uint64_t pos, FSTIter *iter) {
    uint64_t nodeNum = pos >> 8;
    uint8_t cc = pos & 255;

    if (!isTbitSetU(nodeNum, cc)) {
        iter->setVU(level, nodeNum, pos);
        return true;
    }

    level++;
    nodeNum = (iter->positions[level].keyPos < 0) ? childNodeNumU(pos) : (iter->positions[level].keyPos >> 8);

    while (level < cutoff_level_) {
        if (isObitSetU(nodeNum)) {
            iter->setKVU(level, nodeNum, (nodeNum << 8), true);
            return true;
        }
        nextItemU(nodeNum, 0, cc);
        pos = (nodeNum << 8) + cc;
        iter->positions[level].keyPos = pos;

        if (!isTbitSetU(nodeNum, cc)) {
            iter->setVU(level, nodeNum, pos);
            return true;
        }

        level++;
        nodeNum = (iter->positions[level].keyPos < 0) ? childNodeNumU(pos) : (iter->positions[level].keyPos >> 8);
    }

    pos = (iter->positions[level].keyPos < 0) ? childpos(nodeNum) : iter->positions[level].keyPos;
    return nextLeft(level, pos, iter);
}

inline bool FST::nextLeft(int level, uint64_t pos, FSTIter *iter) {
    while (isTbitSet(pos)) {
        iter->positions[level].keyPos = pos;
        level++;
        pos = (iter->positions[level].keyPos < 0) ? childpos(childNodeNum(pos) + childCountU_)
                                                  : iter->positions[level].keyPos;
    }
    iter->setKV(level, pos);
    return true;
}


//******************************************************
// NEXT NODE
//******************************************************
inline bool FST::nextNodeU(int level, uint64_t nodeNum, FSTIter *iter) {
    int cur_level = (level < cutoff_level_) ? level : (cutoff_level_ - 1);
    uint8_t cc = 0;
    uint8_t kc = 0;
    bool inNode = false;

    while (!inNode) {
        nodeNum++;
        if (isObitSetU(nodeNum)) {
            iter->positions[cur_level].keyPos = nodeNum << 8;
            iter->positions[cur_level].isO = true;
        } else {
            nextItemU(nodeNum, 0, cc);
            iter->positions[cur_level].keyPos = (nodeNum << 8) + cc;
        }

        if (cur_level == 0)
            break;

        cur_level--;
        nodeNum = iter->positions[cur_level].keyPos >> 8;
        kc = iter->positions[cur_level].keyPos & 255;

        uint8_t next_kc = (kc == 0 && iter->positions[cur_level].isO) ? kc : (kc + 1);
        iter->positions[cur_level].isO = false;
        inNode = (kc == 255) ? false : nextItemU(nodeNum, next_kc, cc);
    }
    iter->positions[cur_level].keyPos = (nodeNum << 8) + cc;

    while (cur_level < level && cur_level < cutoff_level_) {
        uint64_t pos = iter->positions[cur_level].keyPos;
        nodeNum = pos >> 8;
        cc = pos & 255;
        if (!isTbitSetU(nodeNum, cc)) {
            iter->setVU(cur_level, nodeNum, pos);
            return true;
        }
        cur_level++;
    }

    if (level < cutoff_level_)
        return nextLeftU(level, iter->positions[cur_level].keyPos, iter);
    else
        return false;
}

inline bool FST::nextNode(int level, uint64_t pos, FSTIter *iter) {
    bool inNode = false;
    int cur_level = level - 1;
    while (!inNode && cur_level >= cutoff_level_) {
        iter->positions[cur_level].keyPos++;
        pos = iter->positions[cur_level].keyPos;
        inNode = !isSbitSet(pos);
        cur_level--;
    }

    if (!inNode && cur_level < cutoff_level_) {
        uint64_t nodeNum = iter->positions[cur_level].keyPos >> 8;
        uint8_t kc = iter->positions[cur_level].keyPos & 255;
        uint8_t cc = 0;
        uint8_t next_kc = (kc == 0 && iter->positions[cur_level].isO) ? kc : (kc + 1);
        iter->positions[cur_level].isO = false;

        inNode = (kc == 255) ? false : nextItemU(nodeNum, next_kc, cc);

        if (!inNode) {
            if (nextNodeU(level, (iter->positions[cur_level].keyPos >> 8), iter))
                return true;
        } else {
            iter->positions[cur_level].keyPos = (nodeNum << 8) + cc;
            return nextLeftU(cur_level, iter->positions[cur_level].keyPos, iter);
        }
    }

    cur_level++;
    while (cur_level < level) {
        uint64_t pos = iter->positions[cur_level].keyPos;
        if (!isTbitSet(pos)) {
            iter->setV(cur_level, pos);
            return true;
        }
        cur_level++;
    }

    return nextLeft(level, iter->positions[cur_level].keyPos, iter);
}

//******************************************************
// LOWER BOUND
//******************************************************
bool FST::lowerBound(const uint8_t *key, const int keylen, FSTIter &iter) {
    iter.clear();
    int keypos = 0;
    uint64_t nodeNum = 0;
    uint8_t kc = (uint8_t) key[keypos];
    uint8_t cc = 0;
    uint64_t pos = kc;

    while (keypos < keylen && keypos < cutoff_level_) {
        kc = (uint8_t) key[keypos];
        pos = (nodeNum << 8) + kc;

        __builtin_prefetch(tbitsU_->bits_ + (nodeNum << 2) + (kc >> 6), 0);
        __builtin_prefetch(tbitsU_->rankLUT_ + ((pos + 1) >> 6), 0);

        if (!nextItemU(nodeNum, kc, cc)) { // next char is in next node
            nextNodeU(keypos, nodeNum, &iter);
            return nextLeftU(keypos, iter.positions[keypos].keyPos, &iter);
        }

        if (cc != kc) {
            iter.positions[keypos].keyPos = (nodeNum << 8) + cc;
            return nextLeftU(keypos, iter.positions[keypos].keyPos, &iter);
        }

        iter.positions[keypos].keyPos = pos;

        if (!isTbitSetU(nodeNum, kc)) { // found key terminiation (value)
            iter.len = keypos + 1;
            iter.positions[keypos].valPos = valuePosU(nodeNum, pos);
            return true;
        }

        nodeNum = childNodeNumU(pos);
        keypos++;
    }

    if (keypos < cutoff_level_) {
        pos = nodeNum << 8;
        if (isObitSetU(nodeNum)) {
            iter.setKVU(keypos, nodeNum, pos, true);
            return true;
        }
        keypos--;
        return nextLeftU(keypos, iter.positions[keypos].keyPos, &iter);
    }

    //----------------------------------------------------------
    pos = (cutoff_level_ == 0) ? 0 : childpos(nodeNum);

    bool inNode = true;
    while (keypos < keylen) {
        kc = (uint8_t) key[keypos];

        int nsize = nodeSize(pos);
        inNode = nodeSearch_lowerBound(pos, nsize, kc);

        iter.positions[keypos].keyPos = pos;

        if (!inNode) {
            nextNode(keypos, pos, &iter);
            return nextLeft(keypos, pos, &iter);
        }

        cc = cbytes_[pos];
        if (cc != kc)
            return nextLeft(keypos, pos, &iter);

        if (!isTbitSet(pos)) {
            iter.len = keypos + 1;
            iter.positions[keypos].valPos = valuePos(pos);
            return true;
        }

        pos = childpos(childNodeNum(pos) + childCountU_);
        keypos++;

        __builtin_prefetch(cbytes_ + pos, 0, 1);
        __builtin_prefetch(tbits_->bits_ + (pos >> 6), 0, 1);
        __builtin_prefetch(tbits_->rankLUT_ + ((pos + 1) >> 9), 0);
    }

    if (cbytes_[pos] == TERM && !isTbitSet(pos)) {
        iter.positions[keypos].keyPos = pos;
        iter.len = keypos + 1;
        iter.positions[keypos].valPos = valuePos(pos);
        return true;
    }
    keypos--;
    return nextLeft(keypos, iter.positions[keypos].keyPos, &iter);
}

bool FST::lowerBound(const uint64_t key, FSTIter &iter) {
    uint8_t key_str[8];
    reinterpret_cast<uint64_t *>(key_str)[0] = __builtin_bswap64(key);
    return lowerBound(key_str, 8, iter);
}


//******************************************************
// PRINT
//******************************************************

void FST::print_csv() {
    auto resultU = printU();
    auto resultL = printL();

    std::ofstream out_file;
    std::string out_file_name = "fst.csv";
    out_file.open(out_file_name, std::ios::out | std::ios::app);

    auto upper_nodes_number(0);
    for (int i = 0; i < resultU[2].size(); i++) {
        if (resultU[2][i] == "1") {
            upper_nodes_number++;
        }
    }
    out_file << upper_nodes_number << std::endl;


    for (auto i = 0; i < 4; i++) {
        for (int k = 0; k < resultU[i].size(); k++) {
            out_file << resultU[i][k] << ",";
        }
        for (int k = 0; k < resultL[i].size() - 1; k++) {
            out_file << resultL[i][k] << ",";
        }
        out_file << resultL[i][resultL[i].size() - 1] << std::endl;
    }
    out_file.close();
}

std::vector<std::vector<std::string>> FST::printU() {

    std::vector<std::string> c_bytes;
    std::vector<std::string> t_bits;
    std::vector<std::string> s_bits;
    std::vector<std::string> values;

    for (int i = 0; i < c_lenU_; i += 4) {
        bool first_entry = true;

        for (int j = 0; j < 256; j++) {
            if (isLabelExist(cbitsU_->bits_ + i, (uint8_t) j)) {
                c_bytes.emplace_back(std::to_string((char) j));
                if (first_entry) {
                    s_bits.emplace_back("1");
                    first_entry = false;
                } else {
                    s_bits.emplace_back("0");
                }
            }

        }
    }


    for (int i = 0; i < c_lenU_; i += 4) {
        for (int j = 0; j < 256; j++) {
            if (isLabelExist(tbitsU_->bits_ + i, (uint8_t) j)) {
                t_bits.emplace_back("1");
            } else {
                if (isLabelExist(cbitsU_->bits_ + i, (uint8_t) j)) {
                    t_bits.emplace_back("0");
                }
            }

        }

    }
    assert(t_bits.size() == s_bits.size() && s_bits.size() == c_bytes.size());
    for (int i = 0; i < val_memU_ / 8; i++) {
        values.emplace_back(std::to_string(valuesU_[i]));
    }

    std::vector<std::vector<std::string>> result = {c_bytes, t_bits, s_bits, values};
    return result;
}

std::vector<std::vector<std::string>> FST::printL() {
    std::vector<std::string> c_bytes;
    std::vector<std::string> t_bits;
    std::vector<std::string> s_bits;
    std::vector<std::string> values;

    int c_pos = 0;
    int v_pos = 0;
    int s_pos = 0;
    int value_pos = 0;

    for (int i = 0; i < c_mem_; i++)
        c_bytes.emplace_back(std::to_string(cbytes_[i]));

    for (int i = 0; i < c_mem_; i++) {
        if (readBit(tbits_->getBits()[i / 64], i % 64))
            t_bits.emplace_back("1");
        else
            t_bits.emplace_back("0");
    }


    for (int i = 0; i < c_mem_; i++) {
        if (readBit(sbits_->getBits()[i / 64], i % 64))
            s_bits.emplace_back("1");
        else
            s_bits.emplace_back("0");
    }

    assert(t_bits.size() == s_bits.size() && s_bits.size() == c_bytes.size());

    for (int i = 0; i < val_mem_ / 8; i++) {
        values.emplace_back(std::to_string(values_[i]));
    }
    std::vector<std::vector<std::string>> result = {c_bytes, t_bits, s_bits, values};
    return result;
}


//******************************************************
// ITERATOR
//******************************************************
FSTIter::FSTIter() : index(NULL), len(0), isEnd(false), cBoundU(0), cBound(0), cutoff_level(0), tree_height(0),
                     last_value_pos(0) {}

FSTIter::FSTIter(FST *idx) {
    index = idx;
    tree_height = index->tree_height_;
    cutoff_level = index->cutoff_level_;
    cBoundU = (index->c_lenU_ << 6) - 1;
    cBound = index->cMem() - 1;
    last_value_pos = index->last_value_pos_;

    len = 0;
    isEnd = false;

    for (int i = 0; i < tree_height; i++) {
        Cursor c;
        c.keyPos = -1;
        c.valPos = -1;
        c.isO = false;
        positions.push_back(c);
    }
}

void FSTIter::clear() {
    for (int i = 0; i < tree_height; i++) {
        positions[i].keyPos = -1;
        positions[i].valPos = -1;
        positions[i].isO = false;
    }

    len = 0;
    isEnd = false;
}

inline void FSTIter::setVU(int level, uint64_t nodeNum, uint64_t pos) {
    len = level + 1;
    if (positions[level].valPos < 0)
        positions[level].valPos = index->valuePosU(nodeNum, pos);
    else
        positions[level].valPos++;
}

inline void FSTIter::setKVU(int level, uint64_t nodeNum, uint64_t pos, bool o) {
    positions[level].keyPos = pos;
    positions[level].isO = o;
    len = level + 1;
    if (positions[level].valPos < 0)
        positions[level].valPos = index->valuePosU(nodeNum, pos);
    else
        positions[level].valPos++;
}

inline void FSTIter::setV(int level, uint64_t pos) {
    len = level + 1;
    if (positions[level].valPos < 0)
        positions[level].valPos = index->valuePos(pos);
    else
        positions[level].valPos++;
}

inline void FSTIter::setKV(int level, uint64_t pos) {
    positions[level].keyPos = pos;
    len = level + 1;
    if (positions[level].valPos < 0)
        positions[level].valPos = index->valuePos(pos);
    else
        positions[level].valPos++;
}

//TODO inlining
uint64_t FSTIter::value() {
    if (len <= cutoff_level) {
        __builtin_prefetch(index->valuesU_ + positions[len - 1].valPos + 1, 0, 1);
        return index->valuesU_[positions[len - 1].valPos];
    } else {
        __builtin_prefetch(index->values_ + positions[len - 1].valPos + 1, 0, 1);
        return index->values_[positions[len - 1].valPos];
    }
}

bool FSTIter::operator++(int) {
    if (unlikely(isEnd))
        return false;

    if (unlikely(positions[len - 1].valPos == (0 - last_value_pos) || positions[len - 1].valPos == last_value_pos))
        if ((last_value_pos < 0 && len <= cutoff_level) || (last_value_pos > 0 && len > cutoff_level)) {
            isEnd = true;
            return false;
        }

    uint64_t nodeNum = 0;
    uint8_t kc = 0;
    uint8_t cc = 0;
    bool inNode = true;
    int level = len - 1;
    while (level >= 0) {
        if (level < cutoff_level) {
            if (unlikely(positions[level].keyPos >= cBoundU)) {
                level--;
                continue;
            }

            nodeNum = positions[level].keyPos >> 8;
            kc = positions[level].keyPos & 255;
            uint8_t next_kc = (kc == 0 && positions[level].isO) ? kc : (kc + 1);
            positions[level].isO = false;

            inNode = (kc == 255) ? false : index->nextItemU(nodeNum, next_kc, cc);

            if (!inNode)
                return index->nextNodeU(level, nodeNum, this);

            positions[level].keyPos = (nodeNum << 8) + cc;
            return index->nextLeftU(level, positions[level].keyPos, this);
        } else {
            if (unlikely(positions[level].keyPos >= cBound)) {
                level--;
                continue;
            }

            positions[level].keyPos++;

            if (index->isSbitSet(positions[level].keyPos))
                return index->nextNode(level, positions[level].keyPos, this);

            return index->nextLeft(level, positions[level].keyPos, this);
        }
    }

    return false;
}

//TODO
bool FSTIter::operator--(int) {
    return true;
}

