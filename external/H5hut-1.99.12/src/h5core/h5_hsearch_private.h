#ifndef __H5_HSEARCH_PRIVATE_H
#define __H5_HSEARCH_PRIVATE_H

typedef struct hsearch_data {
	struct _ENTRY *table;
	unsigned int size;
	unsigned int filled;
	int (*compare)(const void*, const void*);
	unsigned int (*compute_hash)(const void*);
	h5_err_t (*free_entry)(const void*);
} h5_hashtable_t; 

/* Action which shall be performed in the call to hsearch.  */
typedef enum {
	H5_FIND,
	H5_ENTER,
	H5_REMOVE
} h5_action_t;

typedef struct h5_entry {
	void* dta;
} h5_entry_t;

/* Reentrant versions which can handle multiple hashing tables at the
   same time.  */
extern h5_err_t
h5priv_hsearch (
	void* item,
	const h5_action_t action,
	void** retval,
	h5_hashtable_t* htab
	);

extern h5_err_t
h5priv_hcreate (
	size_t __nel,
	h5_hashtable_t* __htab,
	int (*compare)(const void*, const void*),
	unsigned int (*compute_hash)(const void*),
	h5_err_t (*free_entry)(const void*)
	);

extern h5_err_t
h5priv_hresize (
	size_t nel,
	h5_hashtable_t* htab
	);

extern h5_err_t
h5priv_hdestroy (
	h5_hashtable_t* __htab
	);

extern h5_err_t
h5priv_hcreate_string_keyed (
	size_t nel,
	h5_hashtable_t* htab,
	h5_err_t (*free_entry)(const void*)
	);

#endif
