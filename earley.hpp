#ifndef EARLEY_HPP
#define EARLEY_HPP

#include <vector>
#include <string>
#include <memory>
#include <map>
#include <set>
#include <deque>
#include <iostream>

/* https://www.cs.unm.edu/~luger/ai-final2/JAVA/CH%2030_The%20Earley%20Parser%20-%20Dynamic%20Programming%20in%20Java.pdf */

#define DOT "•"
#define START "§"

class State;

class RHS
{
public:

    RHS(const std::vector<std::string> &t) : m_terms(t), m_hasDot(false), m_dot(-1)
    {
        for (int i = 0; i < m_terms.size(); i++) {
            if (m_terms[i] == DOT) {
                m_dot = i;
                m_hasDot = true;
                break;
            }
        }
    }

    const int getDotPos() const { return m_dot; }

    const std::vector<std::string> &getTerms() const
    {
        return m_terms;
    }

    const std::string &getPriorToDot() const
    {
        if (m_hasDot && m_dot > 0)
            return m_terms[m_dot-1];

        return m_empty;
    }

    const std::string &getAfterDot() const
    {
        if (m_hasDot && m_dot < m_terms.size() - 1)
            return m_terms[m_dot+1];

        return m_empty;
    }

    bool hasDot() const
    {
        return m_hasDot;
    }

    bool isDotLast() const
    {
        if (m_hasDot)
            return (m_dot == m_terms.size() - 1);

        return false;
    }

    bool isDotFirst() const
    {
        if (m_hasDot)
            return (m_dot == 0);

        return false;
    }

    std::shared_ptr<RHS> addDot() const
    {
        std::vector<std::string> t(m_terms);

        t.insert(t.begin(), DOT);

        return std::make_shared<RHS>(t);
    }

    std::shared_ptr<RHS> addDotLast() const
    {
        std::vector<std::string> t(m_terms);

        t.insert(t.end(), DOT);

        return std::make_shared<RHS>(t);
    }

    std::shared_ptr<RHS> moveDot() const
    {
        std::vector<std::string> t(m_terms);

        for (int i = 0; i < t.size(); i++) {
            if (t[i] == DOT) {
                t[i] = t[i+1];
                t[i+1] = DOT;
                break;
            }
        }

        return std::make_shared<RHS>(t);
    }

    bool operator ==(const RHS& rhs) const
    {
        return m_terms == rhs.m_terms;
    }

private:
    const std::vector<std::string> m_terms;
    bool m_hasDot;
    int m_dot;
    const std::string m_empty = "";
};


class State
{
public:
    typedef std::shared_ptr<State> State_t;

    State(const std::string lhs, std::shared_ptr<RHS> rhs, int i, int j, State_t source=State_t())
        : m_lhs(lhs), m_rhs(rhs), m_i(i), m_j(j)
    {
         if (source)
             m_sources.push_back(source);
    }

    const std::string &getLHS() const { return m_lhs; }

    std::shared_ptr<RHS> getRHS() const { return m_rhs; }

    int getI() const { return m_i; }

    int getJ() const { return m_j; }

    std::string getPriorToDot() const { return m_rhs->getPriorToDot(); }

    const std::string &getAfterDot() const { return m_rhs->getAfterDot(); }

    bool isDotLast() const { return m_rhs->isDotLast(); }

    const std::vector<std::shared_ptr<State> > &getSources() const { return m_sources; }

    void addSources(const std::vector<std::shared_ptr<State> > &sources)
    { 
         m_sources.insert(m_sources.end(), sources.begin(), sources.end()); 
    }

    bool operator ==(const State& rhs) const
    {
        return m_i == rhs.m_i &&
               m_j == rhs.m_j &&
               m_lhs == rhs.m_lhs &&
               *m_rhs == *rhs.m_rhs;
    }

    operator bool() const { return m_lhs.empty(); }

private:
    int m_i, m_j;
    std::string m_lhs;
    std::shared_ptr<RHS> m_rhs;
    std::vector<std::shared_ptr<State> > m_sources;
};

std::ostream& operator <<(std::ostream& os, const State& state)
{
    auto rhs = state.getRHS()->getTerms();
    int rhs_length = rhs.size();
    os << '('  << state.getLHS() << " → ";

    for (int i = 0; i < rhs_length; i++) {
        os << rhs[i];

        if (i < rhs_length - 1) {
            os << ' ';
        }
    }

    os << " [" << state.getI() << " , " << state.getJ() << "]";

    return os << ")";
}


class Chart
{
public:

    Chart() {}

    void addState(std::shared_ptr<State> state)
    {
        bool contains = false;

        for (auto c : m_chart) {
            if (*state == *c) {
                contains = true;
                c->addSources(state->getSources());
                break;
            }
        }

        if (!contains)
            m_chart.push_back(state);

    }

    std::shared_ptr<State> getState(int i) const
    {
        return m_chart.at(i);
    }

    const std::vector<std::shared_ptr<State> > &getStates() const
    {
        return m_chart;
    }

    int size() const
    {
        return m_chart.size();
    }

private:
    std::vector<std::shared_ptr<State> > m_chart;
};

std::ostream& operator <<(std::ostream& os, const Chart &chart)
{
    for (int i = 0; i < chart.size(); i++) {
        os << "S" << i << ": " << *chart.getStates()[i] << std::endl;
    }

    return os;
}


class Grammar
{
public:

    Grammar(const std::string &start) : m_start(start) {}
    virtual ~Grammar() {}

    const std::string &getStart() const { return m_start; }

    std::vector<std::shared_ptr<RHS> > getRHS(const std::string &lhs) const
    {
        std::vector<std::shared_ptr<RHS> > rhs;

        auto range = m_rules.equal_range(lhs);
 
        for (auto i = range.first; i != range.second; i++) {
            rhs.push_back(i->second);
        }

        return rhs;
    }

    void addRule(const std::string &lhs, const std::vector<std::string> &rhs)
    {
        m_rules.insert(std::make_pair(lhs, std::make_shared<RHS>(rhs)));
    }

    void addTerminal(const std::string &t)
    {
        m_terminals.insert(t);
    }

    bool isTerminal(const std::string &token) const
    {
        return m_terminals.find(token) != m_terminals.end();
    }

private:
    std::multimap<const std::string, std::shared_ptr<RHS> > m_rules;
    std::set<std::string> m_terminals;
    std::string m_start;
};


template<class T>
class EarleyParser
{
public:

    EarleyParser(const std::string &start) : m_tokens(nullptr), m_grammar(start) {}

    virtual ~EarleyParser() {}

    Grammar *getGrammar() { return &m_grammar; }

    const std::vector<Chart> *getCharts() const { return &m_charts; }

    std::vector<T> *getTokens() { return m_tokens; }

    /**************************************************************************
     * parse()
     *   This is the main loop for parsing the sentence into the chart. It will
     *   return true if there is at least one successful parse of the sentence.
     *************************************************************************/
    bool parse(std::vector<T> *tokens)
    {
        m_tokens = tokens;

        m_charts = std::vector<Chart>(m_tokens->size() + 1);

        std::vector<std::string> start{DOT, m_grammar.getStart()};

        auto startRHS = std::make_shared<RHS>(start);

        auto startState = std::make_shared<State>(START, startRHS, 0, 0);

        m_charts[0].addState(startState);

        for (int i = 0; i < m_charts.size(); i++) {
            for (int j = 0; j < m_charts[i].size(); j++) {
                auto state = m_charts[i].getState(j);

                auto next_term = state->getAfterDot();

                if (state->isDotLast())
                    completer(state);
                else if (m_grammar.isTerminal(next_term))
                    scanner(state);
                else
                    predictor(state);
            }
        }

        // Determine whether a successful parse.
        std::vector<std::string> fin{m_grammar.getStart(), DOT};
        auto finRHS = std::make_shared<RHS>(fin);
        State finish(START, finRHS, 0, m_tokens->size());

        auto last = m_charts[m_tokens->size()].getState(m_charts[m_tokens->size()].size() - 1);

        return finish == *last;
    }

private:
    /**************************************************************************
     * predictor()
     *   After this function completes all possible states that could 
     *   potentially continue from the state s is added to the charts.
     *************************************************************************/
    void predictor(std::shared_ptr<State> state)
    {
        auto lhs = state->getAfterDot();
        int j = state->getJ();

        for (auto rhs : m_grammar.getRHS(lhs)) {
            m_charts[j].addState(std::make_shared<State>(lhs, rhs->addDot(), j, j, state));
        }
    }

    /**************************************************************************
     * scanner()
     *   After this function completes any rules for the LHS that are 1 term 
     *   only and match the word in the sentence will be added to the chart.
     *************************************************************************/
    void scanner(std::shared_ptr<State> state)
    {
        auto lhs = state->getAfterDot();
        int i = state->getI();
        int j = state->getJ();

        for (auto rhs : m_grammar.getRHS(lhs)) {
            auto terms = rhs->getTerms();

            if (terms.size() == 1 && j < m_tokens->size() && terms[0] == m_tokens->at(j)) {
                m_charts[j+1].addState(std::make_shared<State>(lhs, rhs->addDotLast(), j, j+1, state));
            }
        }
    }

    /**************************************************************************
     * completer()
     *   After this function completes, any state in the i-th chart for which
     *   the string after the dot matches the current state's LHS will be added 
     *   to the j-th chart with the dot moved to the right.
     *************************************************************************/
    void completer(std::shared_ptr<State> state)
    {
        std::string lhs = state->getLHS();

        for (int a = 0; a < m_charts[state->getI()].size(); a++) {
            auto st = m_charts[state->getI()].getState(a);
            auto after = st->getAfterDot();

            if (lhs == after) {
                auto new_state = std::make_shared<State>(st->getLHS(), st->getRHS()->moveDot(), st->getI(), state->getJ(), state);
                m_charts[state->getJ()].addState(new_state);
            }
        }
    }

private:
    Grammar m_grammar;
    std::vector<T> *m_tokens;
    std::vector<Chart> m_charts;
};


class ParseTree
{
public:
    class PTNode;

    typedef std::shared_ptr<ParseTree> ParseTree_t;

    ParseTree() {}

    ParseTree(Grammar *grammar, const std::vector<Chart> *charts) : m_grammar(grammar), m_charts(charts) {}

    ParseTree(const std::string &s)
    {
        m_root = std::make_shared<PTNode>(s, std::shared_ptr<PTNode>());
    }

    ParseTree(const std::string &s, std::shared_ptr<State> st)
    {
        m_root = std::make_shared<PTNode>(s, std::shared_ptr<PTNode>());
        m_stateList.push_back(st);
    }

    ParseTree(std::shared_ptr<PTNode> r) : m_root(r) {}

    std::deque<std::shared_ptr<State> > &stateList() { return m_stateList; }

    std::string prettyPrint() { return m_root->prettyPrint(); }

    std::shared_ptr<PTNode> getRoot() { return m_root; }

    /**************************************************************************
     * getTree()
     *   This starts the parsing of the grammar and charts. It also removes any
     *   duplicate trees that might have been produced during the parsing.
     *   
     *   This is one of the two functions that are public and is the correct way
     *   to get a ParseTree.
    **************************************************************************/
    std::shared_ptr< std::vector<ParseTree_t> > getTree()
    {
        auto lastC = m_charts->at(m_charts->size() - 1);
        auto parse = lastC.getState(lastC.size() - 1);
        std::vector<ParseTree_t> trees;

        std::vector<std::string> fin{m_grammar->getStart(), DOT};
        auto finRHS = std::make_shared<RHS>(fin);
        State finish(START, finRHS, 0, m_charts->size() - 1);

        // If there was a successful parse, find all of the possible parse trees.
        if (*parse == finish) { 
            for (auto state : parse->getSources()) {
                // Find all the trees that could come from this source.
                auto pt = std::make_shared<ParseTree>(parse->getLHS(), parse);

                auto t = parseTree(pt, pt, state);
                trees.insert(trees.end(), t.begin(), t.end());
            }
        }

        // Remove any duplicate trees.
        auto noDups = std::make_shared< std::vector<ParseTree_t> >();

        for (auto pt : trees) {
            bool contains = false;

            for (auto c : *noDups) {
                if (*pt == *c) {
                    contains = true;
                    break;
                }
            }

            if (!contains)
                noDups->push_back(pt);
        }

        return noDups;
    }

private:
    ParseTree_t getParent()
    {
        if (!m_root->getParent())
            return std::make_shared<ParseTree>();

        return std::make_shared<ParseTree>(m_root->getParent());

    }

    ParseTree_t getChild(int i)
    {
        return std::make_shared<ParseTree>(m_root->getChild(i));
    }

    /**************************************************************************
     * addChild()
     *   If the tree does not have a root yet, this child will become the root.
     *   Otherwise it will be added to the root's children.
     *************************************************************************/
    void addChild(const std::string &v)
    {
        if (m_root)
            m_root->addChild(std::make_shared<PTNode>(v, m_root));
        else
            m_root = std::make_shared<PTNode>(v, std::shared_ptr<PTNode>());
    }

    /**************************************************************************
     * copy()
     *   This creates a deep copy of the root and the state list. This is used
     *   when the path we were taking to get to a parse tree has more then one
     *   viable source from our current state.
     *************************************************************************/
    ParseTree_t copy()
    {
        auto t = std::make_shared<ParseTree>(m_root->copy());
        t->stateList().assign(m_stateList.begin(), m_stateList.end());
        return t;
    }

    /**************************************************************************
     * getNodeI()
     *   Do a recursive search to find the node with id i. This assumes that 
     *   there exists a node with such an id.
     *************************************************************************/
    ParseTree_t getNodeI(int i)
    {
        std::shared_ptr<PTNode> p = m_root;

        while (p->id() != i) {
            for (auto c : p->children()) {
                if (i >= c->id()) {
                    p = c;
                    break;
                }
            }
        }

        return std::make_shared<ParseTree>(p);
    }

    /**************************************************************************
     * getRootID()
     *   Returns the id associated with the root.
     *************************************************************************/
    int getRootID()
    {
        return m_root->id();
    }
	
    bool operator ==(const ParseTree &pt)
    {
        return *m_root == *pt.m_root;
    }

    /**************************************************************************
     * parseTree()
     *   Helper function that takes the current ParseTree and the state we are
     *   currently inspecting.
     * 
     *   ParseTree tree is the entire tree that we are working on building up. 
     *   ParseTree child is the current node we are building up that is a 
     *   child of ParseTree tree.
     *   
     *   The tree was started from the right-most leaf. We build up the tree
     *   by adding parents (which becomes the root) and adding all their 
     *   children. 
     *   
     *   An important part of creating the parse tree is determining what is 
     *   the the correct source state. This implementation uses the stateList
     *   for this purpose. The initial state that is in the stateList is
     *   "$ -> S @". 
     *************************************************************************/
    std::vector<ParseTree_t> 
    parseTree(ParseTree_t &tree, const ParseTree_t &child, std::shared_ptr<State> currentState)
    {
        auto trees = std::vector<ParseTree_t>();

        std::vector<std::string> startTerms = {DOT, m_grammar->getStart()};
        auto startRHS = std::make_shared<RHS>(startTerms);
        State start(START, startRHS, 0, 0);

        // If the current state is the start state, we are done. 
        if (*currentState == start) {
            // The only state in the stateList, the currentState and the start state are all same.
            tree->m_stateList.pop_front();
            trees.push_back(tree);
            return trees;
        }

        auto rhs = currentState->getRHS();
        bool addedLHS = false;
        bool alreadyRemoved = false;

        if (rhs->isDotLast()) {
            // This is the currently left-most child of the current node we are working on. 
            child->addChild(currentState->getLHS()); 
            addedLHS = true;
        }

        // If the top state on the stack is the same as the currentState (with 
        //  a dot change) then we remove it. This means that we have backwards
        //  parsed all of the term that are between the dots.
        if (!rhs->isDotLast() && *tree->m_stateList.front()->getRHS() == *currentState->getRHS()->moveDot()) {
            tree->m_stateList.pop_front();
        }

        tree->m_stateList.push_front(currentState);

        auto srcs = currentState->getSources();

        if (m_grammar->isTerminal(currentState->getLHS())) {
            // The currentState is a POS and we don't need these on the stack.
            //  The stack contains states that are not completely backwards 
            //  parsed yet and this one will be complete.
            tree->m_stateList.pop_front();
            auto terms = currentState->getRHS()->getTerms();
            child->getChild(0)->addChild(terms[0]);

            // Select the correct source (the one with the correct i)
            if (srcs.size() == 1) {
                currentState = srcs.at(0);
            } else {
                auto posState = tree->m_stateList.front();
                for (int i = 0; i < srcs.size(); i++){
                    // If the state we are currently looking at matches the 
                    //  the top state in the stateList (with the dot moved 
                    //  one term) then it is the correct source.
                    currentState = srcs.at(i);
                    State moveCurrent(currentState->getLHS(), currentState->getRHS()->moveDot(), 
                                      currentState->getI(), currentState->getJ() + 1);
                    if (moveCurrent == *posState) {
                        break;
                    }
                }
            }

            srcs = currentState->getSources();
            rhs = currentState->getRHS();
            tree->m_stateList.pop_front();

            // If the state is not completely parsed, add it to the stateList.
            if (currentState->getRHS()->getDotPos() > 0)
                tree->m_stateList.push_front(currentState);

            alreadyRemoved = true;
            addedLHS = false;
        }

        // If the dot is first, the state is completely backwards parsed now.
        //  If the state has not already been removed from the stateList, do
        //  so now.
        if (rhs->isDotFirst() && !alreadyRemoved) {
            tree->m_stateList.pop_front();
        }

        // For every source of the updated currentState, we may need to attempt a
        //  backwards parse.
        for (auto nextState: srcs) {
            ParseTree_t treeCopy = tree->copy();
            ParseTree_t childCopy = treeCopy->getNodeI(child->getRootID());
            ParseTree_t nextChild;

            if (rhs->isDotFirst())
                nextChild = childCopy->getParent();
            else if (addedLHS)
                nextChild = childCopy->getChild(0);
            else
                nextChild = childCopy;

            std::string lhs = nextState->getLHS();

            // When the sources list for a state was created during parsing, all of the 
            // potential trees were mixed together. This means that some of the 
            // sources are not valid for the particular parse tree we are creating.
            // The conditions where it could be a source are:
            //    1) The LHS of the potential source is the same as the term prior
            //       to the dot. This source was produced by the predictor step.
            //    2) We have just completed parsing the term prior to the dot. This source
            //       was produced by the completer step. 
            //
            // The scanner step does not need to be handled here due to being handled
            //  when we were handling the POS.
            if (currentState->getRHS()->getPriorToDot() == lhs || 
                 (*tree->m_stateList.front()->getRHS() == *nextState->getRHS()->moveDot() &&
                  tree->m_stateList.front()->getLHS() == nextState->getLHS()))
            {
                auto t = parseTree(treeCopy, nextChild, nextState);
                trees.insert(trees.end(), t.begin(), t.end());
            }
        }

        return trees;
    }

public:
    /**************************************************************************
     * PTNode
     *   This is a helper class. It is the nodes of the parse tree that we 
     *   will be constructing.
     *************************************************************************/
    class PTNode
    {

    public:
        PTNode(const std::string value, std::shared_ptr<PTNode> parent) : m_value(value), m_parent(parent), m_id(ParseTree::m_ID++) {}

        PTNode(const std::string value, std::shared_ptr<PTNode> parent, int id) : m_value(value), m_parent(parent), m_id(id) {}

        int id() { return m_id; }

        const std::string &value() { return m_value; }

        std::shared_ptr<PTNode> getParent() { return m_parent; }

        const std::vector< std::shared_ptr<PTNode> > &children() { return m_children; }

        void addChild(std::shared_ptr<PTNode> c)
        {
            m_children.insert(m_children.begin(), c);
        }

        std::shared_ptr<PTNode> getChild(int i)
        {
            return m_children.at(i);
        }

        void setParent(std::shared_ptr<PTNode> p)
        {
            m_parent = p;
        }

        std::shared_ptr<PTNode> copy()
        {
            auto n = std::make_shared<PTNode>(m_value, m_parent, m_id);

            for (auto child : m_children) {
                auto nc = child->copy();
                n->m_children.push_back(nc);
            }

            return n;
        }

        bool operator ==(PTNode &pt)
        {
            if (m_value != pt.value())
                return false;

            if (m_children.size() != pt.children().size())
                return false;

            for (int i = 0; i < m_children.size(); i++) {
                if (*getChild(i) != *pt.getChild(i))
                    return false;
            }

            return true;
        }

        bool operator !=(PTNode &pt) { return *this == pt ? true : false; }

        std::string prettyPrint(const std::string &offset = "")
        {
            std::string out = offset + m_value + "\n";

            for (auto child : m_children)
                out += child->prettyPrint(offset + "  ");

            return out;
        }

    private:
        int m_id;
        std::string m_value;
        std::shared_ptr<PTNode> m_parent;
        std::vector< std::shared_ptr<PTNode> > m_children;
    };

private:
    static int m_ID;
    std::shared_ptr<PTNode> m_root;
    std::deque<std::shared_ptr<State> > m_stateList;
    Grammar *m_grammar;
    const std::vector<Chart> *m_charts;
};

int ParseTree::m_ID = 0;

#endif // EARLEY_HPP
