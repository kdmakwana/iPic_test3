#include <stdlib.h>
#include <string.h>

#include "h5core/h5_core.h"
#include "h5_core_private.h"


/*
  Functions for ID-list handling.

  ID-lists are data blobs to store ID or indices.
 */
/*
  Allocate new/empty id-list
 */
h5_err_t
h5priv_alloc_idlist (
	h5_loc_idlist_t** list,
	const int32_t	size
	) {
	H5_PRIV_API_ENTER (h5_err_t, "list=%p, size=%d", list, size);
	TRY (*list = h5_calloc (
		     1, sizeof (**list)+size*sizeof ((*list)->items[0])));
	(*list)->size = size;
	H5_PRIV_API_RETURN (H5_SUCCESS);
}

h5_err_t
h5priv_free_idlist (
	h5_loc_idlist_t** list
	) {
	H5_PRIV_API_ENTER (h5_err_t, "list=%p", list);
	if (*list == NULL) H5_PRIV_API_LEAVE (H5_SUCCESS);
	TRY (h5_free (*list));
	*list = NULL;
	H5_PRIV_API_RETURN (H5_SUCCESS);
}

static inline h5_err_t
grow_idlist (
	h5_loc_idlist_t** list,
	size_t new_size
	) {
	H5_PRIV_API_ENTER (h5_err_t, "list=%p, new_size=%zu", list, new_size);
	size_t num_bytes = sizeof (**list) + (new_size-1)*sizeof((*list)->items[0]);
	TRY (*list = h5_alloc (*list, num_bytes));
	(*list)->size = new_size;
	H5_PRIV_FUNC_RETURN (H5_SUCCESS);
}

/*
  Add item to list at position given by \c idx.
*/
h5_loc_idx_t
h5priv_insert_idlist (
	h5_loc_idlist_t** list,
	h5_loc_id_t id,
	h5_loc_idx_t idx
	) {
	H5_PRIV_API_ENTER (h5_err_t,
			   "list=%p, id=%llu, idx=%llu",
			   list,
			   (long long unsigned)id,
			   (long long unsigned)idx);
	if (*list == NULL) {
		TRY (h5priv_alloc_idlist (list, 2));
	} else if ((*list)->num_items == (*list)->size) {
		h5_size_t size = (*list)->size;
		if (size == 0) {
			size = 16;
		} else {
			size *= 2;
		}
		TRY (grow_idlist (list, size));
	}
	h5_loc_idlist_t* l = *list;
	if (idx == -1) {
		idx = l->num_items;
	} else {
		memmove ( 
			&l->items[idx+1],
			&l->items[idx],
			(l->num_items - idx) * sizeof (l->items[0]));
	}
	l->items[idx] = id;
	l->num_items++;
	H5_PRIV_API_RETURN (idx);
}

/*
  Find ID in sorted list
*/
h5_loc_id_t
h5priv_find_idlist (
	h5_loc_idlist_t* list,
	const h5_loc_id_t item
	) {
	H5_PRIV_API_ENTER (h5_err_t,
			   "list=%p, item=%llu",
			   list, (long long unsigned)item);
	if (!list) {
		H5_PRIV_API_LEAVE (-1);
	}
	register size_t low = 0;
	register size_t mid;
	register size_t high = list->num_items - 1;
	register h5_loc_idx_t diff;
	const h5_loc_id_t face_id =  h5tpriv_get_face_id(item);
	const h5_loc_idx_t elem_idx =  h5tpriv_get_elem_idx(item);
	while (low <= high) {
		mid = (low + high) / 2;
		diff = h5tpriv_get_elem_idx(list->items[mid]) - elem_idx;
		// if element indices are equal, we decide on the face indices
		if (diff == 0) {
			diff = h5tpriv_get_face_id (list->items[mid]) - face_id;
		}
           	if ( diff > 0 )
               		high = mid - 1;
           	else if ( diff < 0 )
               		low = mid + 1;
           	else
               		H5_PRIV_API_LEAVE (mid); // found
       	}
       	H5_PRIV_API_RETURN (-(low+1));  // not found
}


/*
  Search in sorted list. If item is not in list, add it.
 */
h5_loc_idx_t
h5priv_search_idlist (
	h5_loc_idlist_t** list,
	h5_loc_id_t item
	) {
	H5_PRIV_API_ENTER (h5_err_t,
			   "list=%p, item=%llu",
			   list, (long long unsigned)item);
	h5_loc_idx_t idx = h5priv_find_idlist (*list, item);
	if (idx < 0) {
		idx = -(idx+1);
		TRY (idx = h5priv_insert_idlist (list, item, idx));
	}
	H5_PRIV_API_RETURN (idx);
}

///////////////////////////////////////////////////////////////////////////////

/*
  Allocate new/empty string-list
 */
h5_err_t
h5priv_alloc_strlist (
	h5_strlist_t** list,
	const h5_size_t	size
	) {
	H5_PRIV_API_ENTER (h5_err_t,
			   "list=%p, size=%llu",
			   list, (long long unsigned)size);
	TRY (*list = h5_calloc (
		     1, sizeof (**list)+size*sizeof ((*list)->items[0])));
	(*list)->size = size;
	H5_PRIV_API_RETURN (H5_SUCCESS);
}

h5_err_t
h5priv_free_strlist (
	h5_strlist_t** list
	) {
	H5_PRIV_API_ENTER (h5_err_t, "list=%p", list);
	if (list == NULL || *list == NULL) H5_PRIV_API_LEAVE (H5_SUCCESS);
	h5_strlist_t* l = *list;
	for (size_t i = 0; i < l->size; i++) {
		TRY (h5_free(l->items[i]));
	} 
	TRY (h5_free (l));
	l = NULL;
	H5_PRIV_API_RETURN (H5_SUCCESS);
}

static inline h5_err_t
grow_strlist (
	h5_strlist_t** list,
	size_t new_size
	) {
	H5_PRIV_API_ENTER (h5_err_t, "list=%p, new_size=%zu", list, new_size);
	size_t num_bytes = sizeof (**list) + (new_size-1)*sizeof((*list)->items[0]);
	TRY (*list = h5_alloc (*list, num_bytes));
	(*list)->size = new_size;
	H5_PRIV_FUNC_RETURN (H5_SUCCESS);
}

/*
  Add item to list at position given by \c idx.
*/
ssize_t
h5priv_insert_strlist (
	h5_strlist_t** list,
	const char const* item,
	size_t idx
	) {
	H5_PRIV_API_ENTER (ssize_t,
			   "list=%p, item=%s, idx=%llu",
			   list, item, (long long unsigned)idx);
	if (*list == NULL) {
		TRY (h5priv_alloc_strlist (list, 2));
	} else if ((*list)->num_items == (*list)->size) {
		h5_size_t size = (*list)->size;
		if (size == 0) {
			size = 16;
		} else {
			size *= 2;
		}
		TRY (grow_strlist (list, size));
	}
	h5_strlist_t* l = *list;
	if (idx == -1) {
		idx = l->num_items;
	} else {
		memmove ( 
			&l->items[idx+1],
			&l->items[idx],
			(l->num_items - idx) * sizeof (l->items[0]));
	}
	TRY (l->items[idx] = h5_strdup(item));
	l->num_items++;
	H5_PRIV_API_RETURN (idx);
}

/*
  Find ID in sorted list
*/
ssize_t
h5priv_find_strlist (
	h5_strlist_t* list,
	const char* const item
	) {
	H5_PRIV_API_ENTER (ssize_t, "list=%p, item=%s", list, item);
	if (!list) {
		H5_PRIV_API_LEAVE (-1);
	}
	register h5_loc_idx_t low = 0;
	register h5_loc_idx_t high = list->num_items - 1;
	register h5_loc_id_t diff;
	register h5_loc_id_t mid;
	while (low <= high) {
		mid = (low + high) / 2;
		diff = strcmp (list->items[mid], item);
           	if ( diff > 0 )
               		high = mid - 1;
           	else if ( diff < 0 )
               		low = mid + 1;
           	else
               		H5_PRIV_API_LEAVE (mid); // found
       	}
       	H5_PRIV_API_RETURN (-(low+1));  // not found
}

/*
  Search in sorted list. If item is not in list, add it.
 */
ssize_t
h5priv_search_strlist (
	h5_strlist_t** list,
	const char* const item
	) {
	H5_PRIV_API_ENTER (h5_err_t, "list=%p, item=%s", list, item);
	h5_loc_idx_t idx = h5priv_find_strlist (*list, item);
	if (idx < 0) {
		idx = -(idx+1);
		TRY (idx = h5priv_insert_strlist (list, item, idx));
	}
	H5_PRIV_API_RETURN (idx);
}

h5_err_t
h5priv_remove_strlist (
	h5_strlist_t* list,
	const char* const item
	) {
	H5_PRIV_API_ENTER (h5_err_t, "list=%p, item=%s", list, item);
	ssize_t idx = h5priv_find_strlist (list, item);
	if (idx < 0) H5_PRIV_API_LEAVE (-1);

	list->num_items--;
	memmove ( 
		&list->items[idx],
		&list->items[idx+1],
		(list->num_items - idx) * sizeof (list->items[0]));
	H5_PRIV_API_RETURN (idx);
}
		

///////////////////////////////////////////////////////////////////////////////

h5_err_t
h5priv_alloc_idxmap (
	h5_idxmap_t* map,
	const h5_size_t	size
	) {
	H5_PRIV_API_ENTER (h5_err_t,
			    "map=%p, size=%llu",
			    map, (long long unsigned)size);
	int new = (map->items == NULL);
	size_t size_in_bytes = size * sizeof (map->items[0]);
	TRY (map->items = h5_alloc (map->items, size_in_bytes));
	map->size = size;
	if (new) map->num_items = 0;
	H5_PRIV_API_RETURN (H5_SUCCESS);
}

h5_err_t
h5priv_insert_idxmap (
	h5_idxmap_t* map,
	h5_glb_idx_t glb_idx,
	h5_loc_idx_t loc_idx
	) {
	H5_PRIV_API_ENTER (h5_err_t,
			   "map=%p, glb_idx=%llu, loc_idx=%llu",
			   map,
			   (long long unsigned)glb_idx,
			   (long long unsigned)loc_idx);
	if (map->num_items == map->size)
		H5_PRIV_API_LEAVE (
			HANDLE_H5_OVERFLOW_ERR (
				(long long)map->size));

	h5_loc_idx_t i = h5priv_search_idxmap (map, glb_idx);
	if (i >= 0)			/* global id already in use ? */
		H5_PRIV_API_LEAVE (-1);

	i = -(i+1);

	memmove ( 
		&map->items[i+1],
		&map->items[i],
		(map->num_items - i) * sizeof(map->items[0]));
	map->items[i].glb_idx = glb_idx;
	map->items[i].loc_idx = loc_idx;
	map->num_items++;

	H5_PRIV_API_RETURN (H5_SUCCESS);
}

/*!

  \ingroup h5_core

  binary search in id map. 

  \return index in array if found, othwise \c -(result+1) is the index
  where \c value must be inserted.

 */
h5_loc_idx_t
h5priv_search_idxmap (
	h5_idxmap_t* map,
	h5_glb_idx_t value
	) {
	H5_PRIV_API_ENTER (h5_err_t,
			   "map=%p, value=%lld",
			   map, (long long)value);
	register h5_loc_idx_t low = 0;
	register h5_loc_idx_t high = map->num_items - 1;
	while (low <= high) {
		register h5_loc_idx_t mid = (low + high) / 2;
		register h5_glb_idx_t diff = map->items[mid].glb_idx - value;
           	if ( diff > 0 )
               		high = mid - 1;
           	else if ( diff < 0 )
               		low = mid + 1;
           	else
               		H5_PRIV_API_LEAVE (mid); // found
       	}
	H5_PRIV_API_RETURN (-(low+1));  // not found
}

static int
cmp_idxmap_items (
	const void* _item1,
	const void* _item2
	) {
	h5_idxmap_el_t* item1 = (h5_idxmap_el_t*)_item1;
	h5_idxmap_el_t* item2 = (h5_idxmap_el_t*)_item2;

	if (item1->glb_idx < item2->glb_idx) return -1;
	if (item1->glb_idx > item2->glb_idx) return 1;
	
	return 0;
}

h5_err_t
h5priv_sort_idxmap (
	h5_idxmap_t* map
	) {
	H5_PRIV_API_ENTER (h5_err_t, "map=%p", map);
	qsort (	map->items, map->num_items, sizeof (map->items[0]),
		cmp_idxmap_items);
	H5_PRIV_API_RETURN (H5_SUCCESS);
}
