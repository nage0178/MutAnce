#include "MutAnce.h"

void list_append(list_t * list, void * data)
{
  if (!list) return;

  /* create list item */
  list_item_t * item = (list_item_t *)xmalloc(sizeof(list_item_t));
  item->data = data;

  /* if list is empty */
  if (!(list->count))
  {
    list->head = list->tail = item;
    list->count = 1;
    item->next = NULL;
    return ;
  }

    list->tail->next = item;
    list->tail = item;
    item->next = NULL;
    list->count++;
}


void list_clear(list_t * list, void (*cb_dealloc)(void *))
{
  list_item_t * head = list->head;

  while (head)
  {
    list_item_t * temp = head;
    head = head->next;
    if (cb_dealloc)
      cb_dealloc(temp->data);
    free(temp);
  }

  list->head = list->tail = NULL;
  list->count = 0;
}

