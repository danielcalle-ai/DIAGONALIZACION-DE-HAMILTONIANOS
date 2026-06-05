#ifndef MALLAS_H
#define MALLAS_H

#include <vector>
#include <memory>
#include <string>

class malla {
public:
    virtual ~malla() = default;
    virtual void construir(int N) = 0;
    virtual const std::vector<double>& getNodos() const = 0;
    virtual const std::vector<double>& getPesos() const = 0;
};

class mallaLegendre : public malla {
private:
    std::vector<double> nodos;
    std::vector<double> pesos;
public:
    void construir(int N) override;
    const std::vector<double>& getNodos() const override;
    const std::vector<double>& getPesos() const override;
};

class mallaLaguerre : public malla {
private:
    std::vector<double> nodos;
    std::vector<double> pesos;
public:
    void construir(int N) override;
    const std::vector<double>& getNodos() const override;
    const std::vector<double>& getPesos() const override;
};

class mallaHermite : public malla {
private:
    std::vector<double> nodos;
    std::vector<double> pesos;
public:
    void construir(int N) override;
    const std::vector<double>& getNodos() const override;
    const std::vector<double>& getPesos() const override;
};

std::unique_ptr<malla> crearmalla(const std::string& tipo_malla);

#endif
