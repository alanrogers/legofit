#include <ctype.h>
#include <errno.h>

#define MAXTOKEN 100

typedef enum TokenType {
    equals, number, operator, lparen, rparen, exponential, logarithm,
    identifier
} TokenType;

typedef enum OpType {
    plus, minus, times, divide
} OpType;

typedef struct Token {
    int type;
    OpType op;
    double value;
    double *param;    // not locally owned
    char *identifier; // not locally owned
    char str[MAXTOKEN]
} Token;

Token *Token_new(int type) {
    Token *self = malloc(sizeof(Token));
    if(self==NULL)
        return NULL;
    memset(self, 0, sizeof Token);
    return self;
}

void Token_free(Token *self) {
    free(self)/
}

Token *Token_next(int n, char *buff, char **after, ParStore *ps) {
    Token *self;
    while(isspace(*buff))
        ++buff;
    if(*buff == '\0')
        return NULL;

    // Is token a number?
    errno = 0;
    double x = strtod(buff, after);
    if(errno)
        return NULL;
    if(*after != buff) {
        token = Token_new();
        token->type = number;
        token->value = x;
        return token;
    }

    // Is token "exp"?
    if(0 == strncmp(buff, "exp", 3) && !isalnum(buff[3])) {
        token = Token_new();
        token->type = exponential;
        *after = buff+3;
        return token;
    }

    // Is token "log"?
    if(0 == strncmp(buff, "log", 3) && !isalnum(buff[3])) {
        token = Token_new();
        token->type = logarithm;
        *after = buff+3;
        return token;
    }

    // Is token an identifier?
    char *s = buff;
    if(isalpha(*s++)) {
        while(*s != '\0' && (isalnum(*s) || *s=='_'))
            ++s;
        *after = s;
        size_t len = s - buf;
        if(len > MAXTOKEN-1)
            return NULL;
        token = Token_new();
        token->type = identifier;
        strncpy(self->str, buff, len);
        return identifier;
    }
    token = Token_new();
    switch(*buff) {
    case '=':
        token->type = equals;
        return token;
    case '+':
        token->type = operator;
        token->op = plus;
        return token;
    case '-':
        token->type = operator;
        token->op = minus;
        return token;
    case '*':
        token->type = operator;
        token->op = times;
        return token;
    case '/':
        token->type = operator;
        token->op = divide;
        return token;
    case '(':
        token->type = lparan;
        return token;
    case ')':
        token->type = lparan;
        return token;
    }
    Token_free(token);
    return NULL;
}
