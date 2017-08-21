// http://www.cs.princeton.edu/~rs/talks/LLRB/LLRB.pdf

typedef enum Color Color;
typedef struct Node Node;
typedef tipId_t key_t;
typedef double val_t;

#define ISRED(node) ((node)!=NULL && (node)->color==red)

enum Color {red, black};

struct Node {
    key_t key;
    val_t value;
    Color color;
    Node *left, *right;
};

Node *Node_new(key_t key, val_t value) {
    Node *self = malloc(sizeof(Node));
    CHECKMEM(self);
    self->key = key;
    self->value = value;
    return self;
}

Node *Node_search(Node *root, key_t key) {
    Node *x = root;
    while(x != NULL) {
        int cmp;
        if(key == x->key)
            return x;
        else if(key < x->key)
            x = x->left;
        else
            x = x->right;
    }
    return NULL;
}

Node *Node_insert(Node *root, key_t key, val_t value) {
    root = Node_insert_r(root, key, value);
    root->color = black;
    return root;
}

Node *Node_insert_r(Node *h, key_t key, val_t value) {
    if(h == NULL)
        return Node_new(key, value);
    if(ISRED(h->left) && ISRED(h->right))
        Node_colorFlip(h);
    if(key == h->key)
        h->value = value;
    else if( key < h->key )
        h->left = Node_insert_r(h->left, key, value);
    else {
        assert(key > h->key);
        h->right = Node_insert_r(h->right, key, value);
    }
    if(ISRED(h->right) == red && !ISRED(h->left))
        h = Node_rotateLeft(h);
    if(ISRED(h->left) && ISRED(h->left->left))
        h = Node_rotateRight(h);
    return h;
}

#if 0
/* The C above was translated from the Java below. */
public class LLRB<Key extends Comparable<Key>, Value>
{
    private static final boolean RED = true;
    private static final boolean BLACK = false;
    private Node root;
    private class Node
    {
        private Key key;
        private Value val;
        private Node left, right;
        private boolean color;
        Node(Key key, Value val)
            {
                this.key = key;
                this.val = val;
                this.color = RED;
            }
    }
    public Value search(Key key)
    {
        Node x = root;
        while (x != null)
            {
                int cmp = key.compareTo(x.key);
                if (cmp == 0) return x.val;
                else if (cmp < 0) x = x.left;
                else if (cmp > 0) x = x.right;
            }
        return null;
    }
    public void insert(Key key, Value value)
    {
        root = insert(root, key, value);
        root.color = BLACK;
    }
    private Node insert(Node h, Key key, Value value)
    {
        if (h == null) return new Node(key, value);
        if (isRed(h.left) && isRed(h.right)) colorFlip(h);
        int cmp = key.compareTo(h.key);
        if (cmp == 0) h.val = value;
        else if (cmp < 0) h.left = insert(h.left, key, value);
        else h.right = insert(h.right, key, value);
        if (isRed(h.right) && !isRed(h.left)) h = rotateLeft(h);
        if (isRed(h.left) && isRed(h.left.left)) h = rotateRight(h);
        return h;
    }
}
#endif
